#include "../utils/pybind_utils.h"
#include "source_base/matrix2d.h"
#include "source_base/math_integral.h"
#include "source_base/math_sphbes.h"
#include "source_base/sph_bessel_tf.h"

#ifdef __MPI
#include "source_base/parallel_common.h"
#endif

#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace pyabacus::utils;

template <typename... Args>
using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

namespace
{
template <typename T>
Matrix2d<T> matrix_from_numpy(
    const py::array_t<T, py::array::c_style | py::array::forcecast>& input_array)
{
    const py::buffer_info input_info = input_array.request();
    if (input_info.ndim != 2)
    {
        throw py::value_error("from_numpy expects a 2-dimensional NumPy array");
    }
    if (input_info.shape[0] <= 0 || input_info.shape[1] <= 0)
    {
        throw py::value_error("from_numpy expects positive matrix dimensions");
    }
    if (input_info.shape[0] > std::numeric_limits<int>::max()
        || input_info.shape[1] > std::numeric_limits<int>::max())
    {
        throw py::value_error("from_numpy matrix dimensions exceed the supported integer range");
    }

    const int number_of_rows = static_cast<int>(input_info.shape[0]);
    const int number_of_cols = static_cast<int>(input_info.shape[1]);
    Matrix2d<T> matrix(number_of_rows, number_of_cols);
    const auto input = input_array.template unchecked<2>();
    for (int row = 0; row < number_of_rows; ++row)
    {
        for (int col = 0; col < number_of_cols; ++col)
        {
            matrix.at_local(row, col) = input(row, col);
        }
    }
    return matrix;
}

template <typename T>
py::array_t<T> matrix_to_numpy(const Matrix2d<T>& matrix)
{
    const Matrix2d<T> global_matrix = matrix.gather();
    const int number_of_rows = global_matrix.get_number_of_global_rows();
    const int number_of_cols = global_matrix.get_number_of_global_cols();
    const std::vector<py::ssize_t> output_shape = {
        static_cast<py::ssize_t>(number_of_rows),
        static_cast<py::ssize_t>(number_of_cols)};
    py::array_t<T> output_array(output_shape);
    auto output = output_array.template mutable_unchecked<2>();
    for (int row = 0; row < number_of_rows; ++row)
    {
        for (int col = 0; col < number_of_cols; ++col)
        {
            output(row, col) = global_matrix.at_local(row, col);
        }
    }
    return output_array;
}

template <typename T>
py::array_t<T> matrix_local_to_numpy(const Matrix2d<T>& matrix)
{
    const int number_of_local_rows = matrix.get_number_of_local_rows();
    const int number_of_local_cols = matrix.get_number_of_local_cols();
    const std::vector<py::ssize_t> output_shape = {
        static_cast<py::ssize_t>(number_of_local_rows),
        static_cast<py::ssize_t>(number_of_local_cols)};
    py::array_t<T> output_array(output_shape);
    auto output = output_array.template mutable_unchecked<2>();
    for (int local_row = 0; local_row < number_of_local_rows; ++local_row)
    {
        for (int local_col = 0; local_col < number_of_local_cols; ++local_col)
        {
            output(local_row, local_col) = matrix.at_local(local_row, local_col);
        }
    }
    return output_array;
}

std::string validate_temporary_numpy_path(const std::string& output_path)
{
    if (output_path.size() < 4 || output_path.compare(output_path.size() - 4, 4, ".npy") != 0)
    {
        throw py::value_error("to_npy output path must end with '.npy'");
    }

    const py::module_ os_module = py::module_::import("os");
    const py::module_ tempfile_module = py::module_::import("tempfile");
    const py::object path_module = os_module.attr("path");
    const std::string resolved_path
        = path_module.attr("realpath")(py::str(output_path)).cast<std::string>();
    const py::object temporary_root_object
        = path_module.attr("realpath")(tempfile_module.attr("gettempdir")());
    const std::string temporary_root = temporary_root_object.cast<std::string>();
    const std::string path_separator = path_module.attr("sep").cast<std::string>();
    std::string temporary_prefix = temporary_root;
    if (!temporary_prefix.empty()
        && temporary_prefix.compare(temporary_prefix.size() - path_separator.size(),
                                    path_separator.size(),
                                    path_separator)
            != 0)
    {
        temporary_prefix += path_separator;
    }
    if (resolved_path == temporary_root
        || resolved_path.compare(0, temporary_prefix.size(), temporary_prefix) != 0)
    {
        throw py::value_error(
            "to_npy output path must be inside the allowed temporary directory");
    }
    return resolved_path;
}

template <typename T>
bool matrix_to_npy(const Matrix2d<T>& matrix,
                   const std::string& output_path,
                   const int root_rank)
{
    static_cast<void>(validate_temporary_numpy_path(output_path));

    int process_rank = 0;
    if (matrix.is_distributed() && !matrix.get_distribution()->is_serial())
    {
#ifdef __MPI
        const MPI_Comm communicator = matrix.get_distribution()->get_mpi_communicator();
        process_rank = Parallel_Common::comm_rank(communicator);
#else
        throw std::logic_error("active distributed matrix requires MPI support");
#endif
    }
    else if (root_rank != 0)
    {
        throw py::value_error("host matrices require root_rank 0");
    }
#ifdef __MPI
    else
    {
        process_rank = Parallel_Common::comm_rank(MPI_COMM_WORLD);
    }
#endif

    if (process_rank != root_rank)
    {
        if (matrix.is_distributed() && !matrix.get_distribution()->is_serial())
        {
            static_cast<void>(matrix.gather_on_root(root_rank));
        }
        return false;
    }

    const Matrix2d<T> global_matrix = matrix.is_distributed()
        ? matrix.gather_on_root(root_rank)
        : matrix;
    py::array_t<T> output_array = matrix_to_numpy(global_matrix);
    py::object save_function = py::module_::import("numpy").attr("save");
    save_function(py::str(output_path),
                  output_array,
                  py::arg("allow_pickle") = false);
    return true;
}

template <typename T>
void bind_matrix2d(py::module& module, const char* class_name)
{
    py::class_<Matrix2d<T>> matrix_class(module, class_name);
    matrix_class
        .def(py::init<const int, const int>(), "rows"_a, "cols"_a)
#ifdef __MPI
        .def(py::init([](const int global_rows,
                         const int global_cols,
                         const int block_size) {
                 return Matrix2d<T>(global_rows,
                                    global_cols,
                                    block_size,
                                    MPI_COMM_WORLD);
             }),
             "global_rows"_a,
             "global_cols"_a,
             "block_size"_a)
        .def(py::init([](const int global_rows,
                         const int global_cols,
                         const int number_of_row_processors,
                         const int number_of_column_processors,
                         const int block_size) {
                 return Matrix2d<T>(global_rows,
                                    global_cols,
                                    number_of_row_processors,
                                    number_of_column_processors,
                                    block_size,
                                    MPI_COMM_WORLD);
             }),
             "global_rows"_a,
             "global_cols"_a,
             "number_of_row_processors"_a,
             "number_of_column_processors"_a,
             "block_size"_a)
#endif
        .def_static("from_numpy", &matrix_from_numpy<T>, "input_array"_a)
        .def("to_numpy", &matrix_to_numpy<T>)
        .def("local_to_numpy", &matrix_local_to_numpy<T>)
        .def("to_npy", &matrix_to_npy<T>, "output_path"_a, "root_rank"_a)
        .def("scatter", [](Matrix2d<T>& matrix, const Matrix2d<T>& global_matrix) {
            matrix.scatter(global_matrix);
        })
        .def("gather", &Matrix2d<T>::gather)
        .def("owns_global_index", &Matrix2d<T>::owns_global_index, "global_row"_a, "global_col"_a)
        .def_property_readonly("global_shape",
                               [](const Matrix2d<T>& matrix) {
                                   return py::make_tuple(
                                       matrix.get_number_of_global_rows(),
                                       matrix.get_number_of_global_cols());
                               })
        .def_property_readonly("local_shape",
                               [](const Matrix2d<T>& matrix) {
                                   return py::make_tuple(
                                       matrix.get_number_of_local_rows(),
                                       matrix.get_number_of_local_cols());
                               })
        .def_property_readonly("shape",
                               [](const Matrix2d<T>& matrix) {
                                   return py::make_tuple(
                                       matrix.get_number_of_global_rows(),
                                       matrix.get_number_of_global_cols());
                               })
        .def_property_readonly("is_distributed", &Matrix2d<T>::is_distributed);
}
} // namespace

void bind_base_math(py::module& m)
{
    // python binding for class Sphbes
    py::class_<ModuleBase::Sphbes>(m, "Sphbes")
        .def(py::init<>())
        .def_static("sphbesj", overload_cast_<const int, const double>()(&ModuleBase::Sphbes::sphbesj), "l"_a, "x"_a)
        .def_static("dsphbesj", overload_cast_<const int, const double>()(&ModuleBase::Sphbes::dsphbesj), "l"_a, "x"_a)
        .def_static("sphbesj",
                    [](const int n, py::array_t<double> r, const double q, const int l, py::array_t<double> jl) {
                        check_1d_array(r, "r");
                        check_1d_array(jl, "jl");
                        ModuleBase::Sphbes::sphbesj(n, get_array_ptr(r), q, l, get_array_ptr(jl));
                    })
        .def_static("dsphbesj",
                    [](const int n, py::array_t<double> r, const double q, const int l, py::array_t<double> djl) {
                        check_1d_array(r, "r");
                        check_1d_array(djl, "djl");
                        ModuleBase::Sphbes::dsphbesj(n, get_array_ptr(r), q, l, get_array_ptr(djl));
                    })
        .def_static("sphbes_zeros", [](const int l, const int n, py::array_t<double> zeros) {
            check_1d_array(zeros, "zeros");
            ModuleBase::Sphbes::sphbes_zeros(l, n, get_array_ptr(zeros));
        });

    // python binding for class Integral
    py::class_<ModuleBase::Integral>(m, "Integral")
        .def(py::init<>())
        .def_static("Simpson_Integral",
                    [](const int mesh, py::array_t<double> func, py::array_t<double> rab, double asum) {
                        check_1d_array(func, "func");
                        check_1d_array(rab, "rab");

                        double isum = asum;
                        ModuleBase::Integral::Simpson_Integral(mesh, get_array_ptr(func), get_array_ptr(rab), isum);
                        return isum;
                    })
        .def_static("Simpson_Integral",
                    [](const int mesh, py::array_t<double> func, const double dr, double asum) {
                        check_1d_array(func, "func");

                        double isum = asum;
                        ModuleBase::Integral::Simpson_Integral(mesh, get_array_ptr(func), dr, isum);
                        return isum;
                    })
        .def_static("Simpson_Integral_0toall",
                    [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum) {
                        check_1d_array(func, "func");
                        check_1d_array(rab, "rab");
                        check_1d_array(asum, "asum");
                        ModuleBase::Integral::Simpson_Integral_0toall(mesh,
                                                                      get_array_ptr(func),
                                                                      get_array_ptr(rab),
                                                                      get_array_ptr(asum));
                    })
        .def_static("Simpson_Integral_alltoinf",
                    [](const int mesh, py::array_t<double> func, py::array_t<double> rab, py::array_t<double> asum) {
                        check_1d_array(func, "func");
                        check_1d_array(rab, "rab");
                        check_1d_array(asum, "asum");
                        ModuleBase::Integral::Simpson_Integral_alltoinf(mesh,
                                                                        get_array_ptr(func),
                                                                        get_array_ptr(rab),
                                                                        get_array_ptr(asum));
                    })
        .def_static("simpson",
                    [](const int n, py::array_t<double> f, const double dx) {
                        check_1d_array(f, "f");
                        return ModuleBase::Integral::simpson(n, get_array_ptr(f), dx);
                    })
        .def_static("simpson",
                    [](const int n, py::array_t<double> f, py::array_t<double> h) {
                        check_1d_array(f, "f");
                        check_1d_array(h, "h");
                        return ModuleBase::Integral::simpson(n, get_array_ptr(f), get_array_ptr(h));
                    })
        .def_static("Gauss_Legendre_grid_and_weight",
                    [](const int n, py::array_t<double> x, py::array_t<double> w) {
                        check_1d_array(x, "x");
                        check_1d_array(w, "w");
                        ModuleBase::Integral::Gauss_Legendre_grid_and_weight(n, get_array_ptr(x), get_array_ptr(w));
                    })
        .def_static(
            "Gauss_Legendre_grid_and_weight",
            [](const double xmin, const double xmax, const int n, py::array_t<double> x, py::array_t<double> w) {
                check_1d_array(x, "x");
                check_1d_array(w, "w");
                ModuleBase::Integral::Gauss_Legendre_grid_and_weight(xmin, xmax, n, get_array_ptr(x), get_array_ptr(w));
            });
    py::class_<ModuleBase::SphericalBesselTransformer>(m, "SphericalBesselTransformer").def(py::init<>());
}

PYBIND11_MODULE(_base_pack, m)
{
    m.doc() = "Submodule for pyabacus: ModuleBase";

    bind_matrix2d<double>(m, "Matrix2dFloat64");
    bind_matrix2d<std::complex<double>>(m, "Matrix2dComplex128");
#ifdef __MPI
    m.attr("has_mpi_support") = true;
#else
    m.attr("has_mpi_support") = false;
#endif

    bind_base_math(m);
}
