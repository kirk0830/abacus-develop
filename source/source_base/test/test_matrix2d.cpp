#ifdef NDEBUG
#undef NDEBUG
#endif

#include "source_base/matrix2d.h"
#include "source_base/module_container/NDArray/distributed_ndarray.h"

#ifdef __MPI
#include "source_base/parallel_common.h"
#endif

#include "gtest/gtest.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

static_assert(std::is_same<Matrix2d<double>, NDArray<double, 2>>::value,
              "Matrix2d<T> must be the same type as NDArray<T, 2>");
static_assert(std::is_same<NDArray<double>, NDArray<double, 0>>::value,
              "NDArray<T> must remain the rank-zero specialization");
static_assert(!std::is_same<NDArray<double>, NDArray<double, 2>>::value,
              "Rank-zero and rank-two NDArray types must be distinct");
static_assert(!std::is_base_of<NDArray<double>, Matrix2d<double>>::value,
              "Matrix2d must not inherit from the rank-zero NDArray");
static_assert(std::is_copy_constructible<Matrix2d<double>>::value,
              "Matrix2d must be copy constructible");
static_assert(std::is_copy_assignable<Matrix2d<double>>::value,
              "Matrix2d must be copy assignable");
static_assert(std::is_nothrow_move_constructible<Matrix2d<double>>::value,
              "Matrix2d move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<Matrix2d<double>>::value,
              "Matrix2d move assignment must be noexcept");

namespace
{
template <typename T>
class HasAt
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().at(0), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasBegin
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().begin(), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasReshape
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().reshape(1, 1), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasShape
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().shape(), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

template <typename T>
class HasSize
{
    template <typename U>
    static auto test(int) -> decltype(std::declval<U&>().size(), std::true_type());

    template <typename>
    static std::false_type test(...);

  public:
    static constexpr bool value = decltype(test<T>(0))::value;
};

static_assert(!HasAt<Matrix2d<double>>::value,
              "Matrix2d must not expose the rank-zero at() interface");
static_assert(!HasBegin<Matrix2d<double>>::value,
              "Matrix2d must not expose rank-zero iterators");
static_assert(!HasReshape<Matrix2d<double>>::value,
              "Matrix2d must not expose rank-zero reshape()");
static_assert(!HasShape<Matrix2d<double>>::value,
              "Matrix2d must not expose a generic shape() container");
static_assert(!HasSize<Matrix2d<double>>::value,
              "Matrix2d must not expose a generic size() container");

void check_host_layout()
{
    Matrix2d<double> matrix(2, 3);

    EXPECT_FALSE(matrix.is_distributed());
    EXPECT_EQ(matrix.get_distribution(), nullptr);
    EXPECT_EQ(matrix.get_number_of_global_rows(), 2);
    EXPECT_EQ(matrix.get_number_of_global_cols(), 3);
    EXPECT_EQ(matrix.get_number_of_local_rows(), 2);
    EXPECT_EQ(matrix.get_number_of_local_cols(), 3);
    EXPECT_EQ(matrix.get_number_of_local_elements(), 6);
    EXPECT_EQ(matrix.get_grid_rows(), 1);
    EXPECT_EQ(matrix.get_grid_cols(), 1);
    EXPECT_EQ(matrix.get_coord_row(), 0);
    EXPECT_EQ(matrix.get_coord_col(), 0);
    EXPECT_EQ(matrix.get_leading_dimension(), 3);
    EXPECT_EQ(matrix.get_descriptor(), nullptr);
    EXPECT_EQ(matrix.get_blacs_context(), -1);
    EXPECT_EQ(matrix.data(), &matrix.at_local(0, 0));

    double value = 1.0;
    for (int row = 0; row < 2; ++row)
    {
        for (int col = 0; col < 3; ++col)
        {
            const std::size_t expected_index = static_cast<std::size_t>(row) * 3U
                + static_cast<std::size_t>(col);
            EXPECT_EQ(&matrix.at_local(row, col), matrix.data() + expected_index);
            EXPECT_EQ(&matrix.at_global(row, col), &matrix.at_local(row, col));
            EXPECT_TRUE(matrix.owns_global_index(row, col));
            matrix.at_local(row, col) = value;
            value += 1.0;
        }
    }

    EXPECT_EQ(matrix.data()[0], 1.0);
    EXPECT_EQ(matrix.data()[1], 2.0);
    EXPECT_EQ(matrix.data()[2], 3.0);
    EXPECT_EQ(matrix.data()[3], 4.0);
    EXPECT_EQ(matrix.data()[4], 5.0);
    EXPECT_EQ(matrix.data()[5], 6.0);

    const Matrix2d<double>& const_matrix = matrix;
    EXPECT_EQ(const_matrix.data(), matrix.data());
    EXPECT_EQ(&const_matrix.at_local(1, 2), matrix.data() + 5);
    EXPECT_EQ(&const_matrix.at_global(1, 2), matrix.data() + 5);
}

void check_distribution_layout(Matrix2d<double>& matrix,
                               const int global_rows,
                               const int global_cols)
{
    const std::shared_ptr<const BlockCyclicDistributionInformation>& distribution
        = matrix.get_distribution();
    ASSERT_NE(distribution, nullptr);
    EXPECT_EQ(matrix.get_number_of_global_rows(), global_rows);
    EXPECT_EQ(matrix.get_number_of_global_cols(), global_cols);
    EXPECT_EQ(matrix.get_number_of_local_rows(), distribution->get_number_of_local_rows());
    EXPECT_EQ(matrix.get_number_of_local_cols(), distribution->get_number_of_local_cols());
    EXPECT_EQ(matrix.get_number_of_local_elements(),
              distribution->get_number_of_local_elements());
    EXPECT_EQ(matrix.get_grid_rows(), distribution->get_grid_rows());
    EXPECT_EQ(matrix.get_grid_cols(), distribution->get_grid_cols());
    EXPECT_EQ(matrix.get_coord_row(), distribution->get_coord_row());
    EXPECT_EQ(matrix.get_coord_col(), distribution->get_coord_col());
    EXPECT_EQ(matrix.get_leading_dimension(), distribution->get_leading_dimension());

    const int local_rows = matrix.get_number_of_local_rows();
    const int local_cols = matrix.get_number_of_local_cols();
    const std::size_t leading_dimension
        = static_cast<std::size_t>(matrix.get_leading_dimension());
    for (int local_row = 0; local_row < local_rows; ++local_row)
    {
        for (int local_col = 0; local_col < local_cols; ++local_col)
        {
            const std::size_t expected_index
                = static_cast<std::size_t>(local_col) * leading_dimension
                + static_cast<std::size_t>(local_row);
            EXPECT_EQ(&matrix.at_local(local_row, local_col), matrix.data() + expected_index);
        }
    }

    for (int global_row = 0; global_row < global_rows; ++global_row)
    {
        for (int global_col = 0; global_col < global_cols; ++global_col)
        {
            const bool expected_owner = distribution->owns_index(global_row, global_col);
            EXPECT_EQ(matrix.owns_global_index(global_row, global_col), expected_owner);
            if (!expected_owner)
            {
                continue;
            }

            const int local_row = distribution->global_to_local_row(global_row);
            const int local_col = distribution->global_to_local_col(global_col);
            ASSERT_NE(local_row, -1);
            ASSERT_NE(local_col, -1);
            EXPECT_EQ(&matrix.at_global(global_row, global_col),
                      &matrix.at_local(local_row, local_col));
            const double value = static_cast<double>(global_row * global_cols + global_col);
            matrix.at_global(global_row, global_col) = value;
            EXPECT_EQ(matrix.at_local(local_row, local_col), value);
        }
    }

    const Matrix2d<double>& const_matrix = matrix;
    EXPECT_EQ(const_matrix.data(), matrix.data());
}

#ifdef __MPI
void check_mpi_descriptor(Matrix2d<double>& matrix,
                          const int global_rows,
                          const int global_cols,
                          const int block_size)
{
    constexpr int kDescriptorMatrixTypeIndex = 0;
    constexpr int kDescriptorBlacsContextIndex = 1;
    constexpr int kDescriptorGlobalRowsIndex = 2;
    constexpr int kDescriptorGlobalColsIndex = 3;
    constexpr int kDescriptorBlockRowsIndex = 4;
    constexpr int kDescriptorBlockColsIndex = 5;
    constexpr int kDescriptorLeadingDimIndex = 8;

    const int* const descriptor = matrix.get_descriptor();
    ASSERT_NE(descriptor, nullptr);
    EXPECT_EQ(descriptor[kDescriptorMatrixTypeIndex], 1);
    EXPECT_EQ(descriptor[kDescriptorBlacsContextIndex], matrix.get_blacs_context());
    EXPECT_EQ(descriptor[kDescriptorGlobalRowsIndex], global_rows);
    EXPECT_EQ(descriptor[kDescriptorGlobalColsIndex], global_cols);
    EXPECT_EQ(descriptor[kDescriptorBlockRowsIndex], block_size);
    EXPECT_EQ(descriptor[kDescriptorBlockColsIndex], block_size);
    EXPECT_EQ(descriptor[kDescriptorLeadingDimIndex], matrix.get_leading_dimension());
    EXPECT_GE(matrix.get_blacs_context(), 0);
    EXPECT_EQ(matrix.get_grid_rows() * matrix.get_grid_cols(),
              Parallel_Common::comm_size(MPI_COMM_WORLD));
}
#endif
} // namespace

TEST(Matrix2dTest, HostConstructionAndRowMajorLayout)
{
    check_host_layout();
}

TEST(Matrix2dTest, HostCopyDeepCopiesStorage)
{
    Matrix2d<double> source(2, 2);
    source.at_local(0, 0) = 1.0;
    source.at_local(0, 1) = 2.0;
    source.at_local(1, 0) = 3.0;
    source.at_local(1, 1) = 4.0;

    Matrix2d<double> copy(source);
    EXPECT_EQ(copy.get_distribution(), nullptr);
    EXPECT_NE(copy.data(), source.data());
    EXPECT_EQ(copy.at_local(1, 1), 4.0);
    copy.at_local(1, 1) = 9.0;
    EXPECT_EQ(source.at_local(1, 1), 4.0);

    Matrix2d<double> assigned(1, 1);
    assigned = source;
    EXPECT_NE(assigned.data(), source.data());
    EXPECT_EQ(assigned.at_local(0, 1), 2.0);
}

TEST(Matrix2dTest, SerialDistributionHandleUsesColumnMajorLayout)
{
    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = make_frozen_distribution(3, 2);
    Matrix2d<double> matrix(distribution);

    EXPECT_TRUE(matrix.is_distributed());
    EXPECT_EQ(matrix.get_distribution().get(), distribution.get());
    EXPECT_EQ(matrix.get_blacs_context(), -1);
    check_distribution_layout(matrix, 3, 2);
}

TEST(Matrix2dTest, CopiesShareDistributionAndDeepCopyStorage)
{
    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = make_frozen_distribution(3, 2);
    Matrix2d<double> source(distribution);
    source.at_local(0, 0) = 7.0;

    Matrix2d<double> copy(source);
    EXPECT_EQ(copy.get_distribution().get(), distribution.get());
    EXPECT_EQ(source.get_distribution().get(), copy.get_distribution().get());
    EXPECT_NE(copy.data(), source.data());
    EXPECT_EQ(copy.at_local(0, 0), 7.0);
    copy.at_local(0, 0) = 8.0;
    EXPECT_EQ(source.at_local(0, 0), 7.0);
}

TEST(Matrix2dTest, MatricesFromSameHandleShareContext)
{
    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = make_frozen_distribution(5, 4);
    Matrix2d<double> first(distribution);
    Matrix2d<double> second(distribution);

    EXPECT_EQ(first.get_distribution().get(), second.get_distribution().get());
    EXPECT_EQ(first.get_blacs_context(), second.get_blacs_context());
    EXPECT_NE(first.data(), second.data());
}

TEST(Matrix2dTest, RejectsNullDistribution)
{
    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution;
    EXPECT_THROW(
        {
            Matrix2d<double> matrix(distribution);
            static_cast<void>(matrix);
        },
        std::invalid_argument);
}

#ifdef __MPI
TEST(Matrix2dTest, BlockSizeConstructorAtLocalAndAtGlobal)
{
    const int global_rows = 8;
    const int global_cols = 6;
    const int block_size = 2;
    Matrix2d<double> matrix(global_rows, global_cols, block_size, MPI_COMM_WORLD);

    check_mpi_descriptor(matrix, global_rows, global_cols, block_size);
    check_distribution_layout(matrix, global_rows, global_cols);
}

TEST(Matrix2dTest, UnitBlockConstructorAtLocalAndAtGlobal)
{
    const int global_rows = 5;
    const int global_cols = 3;
    Matrix2d<double> matrix(global_rows, global_cols, MPI_COMM_WORLD);

    check_mpi_descriptor(matrix, global_rows, global_cols, 1);
    check_distribution_layout(matrix, global_rows, global_cols);
}

TEST(Matrix2dTest, FrozenHandleConstructorSharesBlacsContext)
{
    const int global_rows = 8;
    const int global_cols = 6;
    const int block_size = 2;
    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = make_frozen_distribution(global_rows, global_cols, block_size, MPI_COMM_WORLD);
    Matrix2d<double> first(distribution);
    Matrix2d<double> second(distribution);

    const int distribution_context = distribution->get_blacs_context();
    EXPECT_GE(distribution_context, 0);
    EXPECT_EQ(first.get_distribution().get(), distribution.get());
    EXPECT_EQ(second.get_distribution().get(), distribution.get());
    EXPECT_EQ(first.get_blacs_context(), distribution_context);
    EXPECT_EQ(second.get_blacs_context(), distribution_context);
    check_mpi_descriptor(first, global_rows, global_cols, block_size);
}
#endif

int main(int argc, char** argv)
{
#ifdef __MPI
    MPI_Init(&argc, &argv);
#endif
    testing::InitGoogleTest(&argc, argv);
    const int result = RUN_ALL_TESTS();
#ifdef __MPI
    MPI_Finalize();
#endif
    return result;
}
