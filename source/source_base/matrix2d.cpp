#include "source_base/matrix2d.h"

#include "source_base/parallel_common.h"
#ifdef __MPI
#include "source_base/module_external/scalapack_connector.h"
#endif

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{
constexpr int kLocalMetadataSize = 6;

template <typename T>
T conjugate_element(const T& value)
{
    return value;
}

template <>
std::complex<float> conjugate_element<std::complex<float>>(
    const std::complex<float>& value)
{
    return std::conj(value);
}

template <>
std::complex<double> conjugate_element<std::complex<double>>(
    const std::complex<double>& value)
{
    return std::conj(value);
}

template <typename T>
double real_element(const T& value)
{
    return static_cast<double>(value);
}

template <>
double real_element<std::complex<float>>(const std::complex<float>& value)
{
    return static_cast<double>(value.real());
}

template <>
double real_element<std::complex<double>>(const std::complex<double>& value)
{
    return value.real();
}

template <typename T>
double imaginary_element(const T&)
{
    return 0.0;
}

template <>
double imaginary_element<std::complex<float>>(const std::complex<float>& value)
{
    return static_cast<double>(value.imag());
}

template <>
double imaginary_element<std::complex<double>>(const std::complex<double>& value)
{
    return value.imag();
}

template <typename T>
double squared_element(const T& value)
{
    const double real_value = static_cast<double>(value);
    return real_value * real_value;
}

template <>
double squared_element<std::complex<float>>(const std::complex<float>& value)
{
    const double real_value = static_cast<double>(value.real());
    const double imaginary_value = static_cast<double>(value.imag());
    return real_value * real_value + imaginary_value * imaginary_value;
}

template <>
double squared_element<std::complex<double>>(const std::complex<double>& value)
{
    return value.real() * value.real() + value.imag() * value.imag();
}

template <typename T>
Matrix2d<double> extract_real_or_imaginary(const Matrix2d<T>& matrix,
                                           const bool extract_real)
{
    Matrix2d<double> result = matrix.is_distributed()
        ? Matrix2d<double>(matrix.get_distribution())
        : Matrix2d<double>(matrix.get_number_of_global_rows(),
                           matrix.get_number_of_global_cols());

    if (!matrix.is_distributed())
    {
        for (int row = 0; row < matrix.get_number_of_global_rows(); ++row)
        {
            for (int col = 0; col < matrix.get_number_of_global_cols(); ++col)
            {
                const T& value = matrix.at_local(row, col);
                result.at_local(row, col)
                    = extract_real ? real_element(value) : imaginary_element(value);
            }
        }
        return result;
    }

    for (int local_row = 0; local_row < matrix.get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < matrix.get_number_of_local_cols(); ++local_col)
        {
            const T& value = matrix.at_local(local_row, local_col);
            result.at_local(local_row, local_col)
                = extract_real ? real_element(value) : imaginary_element(value);
        }
    }
    return result;
}

template <typename T>
void write_real_or_imaginary(const Matrix2d<T>& matrix,
                             Matrix2d<double>& output,
                             const bool extract_real)
{
    if (!matrix.is_distributed())
    {
        for (int row = 0; row < matrix.get_number_of_global_rows(); ++row)
        {
            for (int col = 0; col < matrix.get_number_of_global_cols(); ++col)
            {
                const T& value = matrix.at_local(row, col);
                output.at_local(row, col)
                    = extract_real ? real_element(value) : imaginary_element(value);
            }
        }
        return;
    }

    for (int local_row = 0; local_row < matrix.get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < matrix.get_number_of_local_cols(); ++local_col)
        {
            const T& value = matrix.at_local(local_row, local_col);
            output.at_local(local_row, local_col)
                = extract_real ? real_element(value) : imaginary_element(value);
        }
    }
}

#ifdef __MPI
int checked_mpi_count(const std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        throw std::invalid_argument("Matrix2d collective payload exceeds MPI count range");
    }
    return static_cast<int>(count);
}

std::vector<int> make_byte_displacements(const std::vector<int>& counts)
{
    std::vector<int> displacements(counts.size(), 0);
    int displacement = 0;
    for (std::size_t rank = 0; rank < counts.size(); ++rank)
    {
        displacements[rank] = displacement;
        if (counts[rank] > std::numeric_limits<int>::max() - displacement)
        {
            throw std::invalid_argument("Matrix2d collective payload size overflow");
        }
        displacement += counts[rank];
    }
    return displacements;
}

std::vector<int> allgather_byte_counts(const int local_count, const MPI_Comm communicator)
{
    const int number_of_ranks = Parallel_Common::comm_size(communicator);
    std::vector<int> receive_counts(number_of_ranks, checked_mpi_count(sizeof(int)));
    const std::vector<int> receive_displacements
        = make_byte_displacements(receive_counts);
    std::vector<int> counts(number_of_ranks, 0);

    Parallel_Common::allgatherv_bytes(&local_count,
                                      checked_mpi_count(sizeof(local_count)),
                                      counts.data(),
                                      receive_counts.data(),
                                      receive_displacements.data(),
                                      communicator);
    return counts;
}

template <typename T>
std::vector<int> make_local_metadata(const Matrix2d<T>& matrix)
{
    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = matrix.get_distribution();
    std::vector<int> metadata;
    metadata.reserve(kLocalMetadataSize + matrix.get_number_of_local_rows()
                     + matrix.get_number_of_local_cols());
    metadata.push_back(matrix.get_number_of_local_rows());
    metadata.push_back(matrix.get_number_of_local_cols());
    metadata.push_back(matrix.get_number_of_global_rows());
    metadata.push_back(matrix.get_number_of_global_cols());
    metadata.push_back(distribution->get_coord_row());
    metadata.push_back(distribution->get_coord_col());
    for (int local_row = 0; local_row < matrix.get_number_of_local_rows(); ++local_row)
    {
        metadata.push_back(distribution->local_to_global_row(local_row));
    }
    for (int local_col = 0; local_col < matrix.get_number_of_local_cols(); ++local_col)
    {
        metadata.push_back(distribution->local_to_global_col(local_col));
    }
    return metadata;
}

template <typename T>
std::vector<unsigned char> make_local_payload(const Matrix2d<T>& matrix)
{
    const std::vector<int> metadata = make_local_metadata(matrix);
    const std::size_t metadata_bytes = metadata.size() * sizeof(int);
    const std::size_t data_elements = static_cast<std::size_t>(matrix.get_number_of_local_rows())
        * static_cast<std::size_t>(matrix.get_number_of_local_cols());
    const std::size_t data_bytes = data_elements * sizeof(T);
    std::vector<unsigned char> payload(metadata_bytes + data_bytes);

    std::memcpy(payload.data(), metadata.data(), metadata_bytes);
    for (int local_row = 0; local_row < matrix.get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < matrix.get_number_of_local_cols(); ++local_col)
        {
            const std::size_t local_row_index = static_cast<std::size_t>(local_row);
            const std::size_t local_col_index = static_cast<std::size_t>(local_col);
            const std::size_t element_index
                = local_row_index * static_cast<std::size_t>(matrix.get_number_of_local_cols())
                + local_col_index;
            const T& value = matrix.at_local(local_row, local_col);
            std::memcpy(payload.data() + metadata_bytes + element_index * sizeof(T),
                        &value,
                        sizeof(T));
        }
    }
    return payload;
}

struct LocalPayloadView
{
    int local_rows = 0;
    int local_cols = 0;
    int global_rows = 0;
    int global_cols = 0;
    int coord_row = 0;
    int coord_col = 0;
    std::vector<int> row_axes;
    std::vector<int> col_axes;
    const unsigned char* data = nullptr;
};

template <typename T>
LocalPayloadView view_local_payload(const std::vector<unsigned char>& payload)
{
    if (payload.size() < kLocalMetadataSize * sizeof(int))
    {
        throw std::invalid_argument("Matrix2d collective payload is truncated");
    }

    int header[kLocalMetadataSize] = {};
    std::memcpy(header, payload.data(), sizeof(header));
    LocalPayloadView view;
    view.local_rows = header[0];
    view.local_cols = header[1];
    view.global_rows = header[2];
    view.global_cols = header[3];
    view.coord_row = header[4];
    view.coord_col = header[5];

    if (view.local_rows < 0 || view.local_cols < 0 || view.global_rows <= 0
        || view.global_cols <= 0)
    {
        throw std::invalid_argument("Matrix2d collective payload has invalid dimensions");
    }

    const std::size_t metadata_count
        = static_cast<std::size_t>(kLocalMetadataSize) + view.local_rows + view.local_cols;
    const std::size_t metadata_bytes = metadata_count * sizeof(int);
    const std::size_t data_elements
        = static_cast<std::size_t>(view.local_rows) * static_cast<std::size_t>(view.local_cols);
    const std::size_t expected_bytes = metadata_bytes + data_elements * sizeof(T);
    if (payload.size() != expected_bytes)
    {
        throw std::invalid_argument("Matrix2d collective payload has an invalid size");
    }

    std::vector<int> metadata(metadata_count);
    std::memcpy(metadata.data(), payload.data(), metadata_bytes);
    view.row_axes.assign(metadata.begin() + kLocalMetadataSize,
                         metadata.begin() + kLocalMetadataSize + view.local_rows);
    view.col_axes.assign(metadata.begin() + kLocalMetadataSize + view.local_rows,
                         metadata.end());
    view.data = payload.data() + metadata_bytes;
    return view;
}

template <typename T>
Matrix2d<T> rebuild_global_matrix(
    const std::vector<std::vector<unsigned char>>& payloads)
{
    std::vector<LocalPayloadView> views;
    views.reserve(payloads.size());
    int global_rows = 0;
    int global_cols = 0;
    for (std::size_t rank = 0; rank < payloads.size(); ++rank)
    {
        views.push_back(view_local_payload<T>(payloads[rank]));
        const LocalPayloadView& view = views.back();
        if (rank == 0)
        {
            global_rows = view.global_rows;
            global_cols = view.global_cols;
        }
        else if (view.global_rows != global_rows || view.global_cols != global_cols)
        {
            throw std::invalid_argument("Matrix2d ranks disagree on global shape");
        }
    }

    Matrix2d<T> result(global_rows, global_cols);
    for (std::size_t rank = 0; rank < views.size(); ++rank)
    {
        const LocalPayloadView& view = views[rank];
        for (int local_row = 0; local_row < view.local_rows; ++local_row)
        {
            const int global_row = view.row_axes[local_row];
            if (global_row < 0 || global_row >= global_rows)
            {
                throw std::invalid_argument("Matrix2d row axis is outside global shape");
            }
            for (int local_col = 0; local_col < view.local_cols; ++local_col)
            {
                const int global_col = view.col_axes[local_col];
                if (global_col < 0 || global_col >= global_cols)
                {
                    throw std::invalid_argument("Matrix2d column axis is outside global shape");
                }
                const std::size_t element_index
                    = static_cast<std::size_t>(local_row)
                        * static_cast<std::size_t>(view.local_cols)
                    + static_cast<std::size_t>(local_col);
                T value;
                std::memcpy(&value, view.data + element_index * sizeof(T), sizeof(T));
                result.at_local(global_row, global_col) = value;
            }
        }
    }
    return result;
}

template <typename T>
std::vector<std::vector<unsigned char>> allgather_payloads(const Matrix2d<T>& matrix,
                                                          const MPI_Comm communicator)
{
    const std::vector<unsigned char> local_payload = make_local_payload(matrix);
    const int local_count = checked_mpi_count(local_payload.size());
    std::vector<int> payload_counts = allgather_byte_counts(local_count, communicator);
    const std::vector<int> payload_displacements
        = make_byte_displacements(payload_counts);
    const int total_bytes = payload_displacements.back() + payload_counts.back();
    std::vector<unsigned char> received(total_bytes);

    Parallel_Common::allgatherv_bytes(local_payload.data(),
                                      local_count,
                                      received.data(),
                                      payload_counts.data(),
                                      payload_displacements.data(),
                                      communicator);

    std::vector<std::vector<unsigned char>> payloads(payload_counts.size());
    for (std::size_t rank = 0; rank < payload_counts.size(); ++rank)
    {
        const std::size_t begin = static_cast<std::size_t>(payload_displacements[rank]);
        const std::size_t end
            = begin + static_cast<std::size_t>(payload_counts[rank]);
        payloads[rank].assign(received.begin() + begin, received.begin() + end);
    }
    return payloads;
}

template <typename T>
std::unique_ptr<Matrix2d<T>> gather_payloads_on_root(const Matrix2d<T>& matrix,
                                                     const int root_rank,
                                                     const MPI_Comm communicator)
{
    const int number_of_ranks = Parallel_Common::comm_size(communicator);
    const int rank = Parallel_Common::comm_rank(communicator);
    if (root_rank < 0 || root_rank >= number_of_ranks)
    {
        throw std::invalid_argument("gather root rank is outside the communicator");
    }

    const std::vector<unsigned char> local_payload = make_local_payload(matrix);
    const int local_count = checked_mpi_count(local_payload.size());
    std::vector<int> payload_counts = allgather_byte_counts(local_count, communicator);
    const std::vector<int> payload_displacements
        = make_byte_displacements(payload_counts);
    const int total_bytes = payload_displacements.back() + payload_counts.back();
    std::vector<unsigned char> received;
    if (rank == root_rank)
    {
        received.resize(total_bytes);
    }

    unsigned char* const receive_buffer
        = rank == root_rank ? received.data() : nullptr;
    Parallel_Common::gatherv_bytes(local_payload.data(),
                                   local_count,
                                   receive_buffer,
                                   payload_counts.data(),
                                   payload_displacements.data(),
                                   root_rank,
                                   communicator);
    if (rank != root_rank)
    {
        return std::unique_ptr<Matrix2d<T>>();
    }

    std::vector<std::vector<unsigned char>> payloads(payload_counts.size());
    for (std::size_t source_rank = 0; source_rank < payload_counts.size(); ++source_rank)
    {
        const std::size_t begin
            = static_cast<std::size_t>(payload_displacements[source_rank]);
        const std::size_t end
            = begin + static_cast<std::size_t>(payload_counts[source_rank]);
        payloads[source_rank].assign(received.begin() + begin, received.begin() + end);
    }
    return std::unique_ptr<Matrix2d<T>>(new Matrix2d<T>(rebuild_global_matrix<T>(payloads)));
}

template <typename T>
void scatter_from_host(const Matrix2d<T>& global_matrix, Matrix2d<T>& local_matrix)
{
    if (global_matrix.is_distributed())
    {
        throw std::invalid_argument("scatter requires a host global matrix");
    }

    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = local_matrix.get_distribution();
    const MPI_Comm communicator = distribution->get_mpi_communicator();
    if (communicator == MPI_COMM_NULL)
    {
        throw std::invalid_argument("scatter requires a valid MPI communicator");
    }

    const std::vector<int> local_meta = make_local_metadata(local_matrix);
    const int local_metadata_bytes
        = checked_mpi_count(local_meta.size() * sizeof(int));
    std::vector<int> metadata_counts
        = allgather_byte_counts(local_metadata_bytes, communicator);
    const std::vector<int> metadata_displacements
        = make_byte_displacements(metadata_counts);
    const int metadata_total_bytes
        = metadata_displacements.back() + metadata_counts.back();
    if (metadata_total_bytes % checked_mpi_count(sizeof(int)) != 0)
    {
        throw std::invalid_argument("Matrix2d scatter metadata has an invalid byte count");
    }
    std::vector<int> all_metadata(metadata_total_bytes / sizeof(int));

    Parallel_Common::allgatherv_bytes(local_meta.data(),
                                      local_metadata_bytes,
                                      all_metadata.data(),
                                      metadata_counts.data(),
                                      metadata_displacements.data(),
                                      communicator);

    const int number_of_ranks = Parallel_Common::comm_size(communicator);
    const int rank = Parallel_Common::comm_rank(communicator);
    std::vector<int> send_counts(number_of_ranks, 0);
    std::vector<LocalPayloadView> views;
    views.reserve(number_of_ranks);
    std::size_t metadata_offset = 0;
    const int global_rows = local_matrix.get_number_of_global_rows();
    const int global_cols = local_matrix.get_number_of_global_cols();
    for (int source_rank = 0; source_rank < number_of_ranks; ++source_rank)
    {
        if (metadata_offset + kLocalMetadataSize > all_metadata.size())
        {
            throw std::invalid_argument("Matrix2d scatter metadata is truncated");
        }
        LocalPayloadView view;
        view.local_rows = all_metadata[metadata_offset];
        view.local_cols = all_metadata[metadata_offset + 1];
        view.global_rows = all_metadata[metadata_offset + 2];
        view.global_cols = all_metadata[metadata_offset + 3];
        view.coord_row = all_metadata[metadata_offset + 4];
        view.coord_col = all_metadata[metadata_offset + 5];
        metadata_offset += kLocalMetadataSize;
        const std::size_t axes_count
            = static_cast<std::size_t>(view.local_rows + view.local_cols);
        if (metadata_offset + axes_count > all_metadata.size())
        {
            throw std::invalid_argument("Matrix2d scatter axes are truncated");
        }
        view.row_axes.assign(all_metadata.begin() + metadata_offset,
                             all_metadata.begin() + metadata_offset + view.local_rows);
        view.col_axes.assign(all_metadata.begin() + metadata_offset + view.local_rows,
                             all_metadata.begin() + metadata_offset + view.local_rows
                                 + view.local_cols);
        metadata_offset += axes_count;

        if (view.local_rows < 0 || view.local_cols < 0
            || view.global_rows != global_rows || view.global_cols != global_cols)
        {
            throw std::invalid_argument("Matrix2d ranks disagree on scatter global shape");
        }
        const std::size_t local_elements
            = static_cast<std::size_t>(view.local_rows)
            * static_cast<std::size_t>(view.local_cols);
        send_counts[source_rank] = checked_mpi_count(local_elements * sizeof(T));
        views.push_back(view);
    }
    if (metadata_offset != all_metadata.size())
    {
        throw std::invalid_argument("Matrix2d scatter metadata has trailing entries");
    }
    if (global_matrix.get_number_of_global_rows() != global_rows
        || global_matrix.get_number_of_global_cols() != global_cols)
    {
        throw std::invalid_argument("scatter global matrix shape does not match distribution");
    }

    const std::vector<int> send_displacements = make_byte_displacements(send_counts);
    std::vector<unsigned char> send_buffer;
    if (rank == 0)
    {
        const int total_bytes = send_displacements.back() + send_counts.back();
        send_buffer.resize(total_bytes);
        for (int destination_rank = 0; destination_rank < number_of_ranks;
             ++destination_rank)
        {
            const LocalPayloadView& view = views[destination_rank];
            for (int local_row = 0; local_row < view.local_rows; ++local_row)
            {
                const int global_row = view.row_axes[local_row];
                for (int local_col = 0; local_col < view.local_cols; ++local_col)
                {
                    const int global_col = view.col_axes[local_col];
                    const std::size_t element_index
                        = static_cast<std::size_t>(local_row)
                            * static_cast<std::size_t>(view.local_cols)
                        + static_cast<std::size_t>(local_col);
                    const T& value = global_matrix.at_local(global_row, global_col);
                    std::memcpy(send_buffer.data()
                                    + send_displacements[destination_rank]
                                    + element_index * sizeof(T),
                                &value,
                                sizeof(T));
                }
            }
        }
    }

    const std::size_t receive_elements
        = static_cast<std::size_t>(local_matrix.get_number_of_local_rows())
        * static_cast<std::size_t>(local_matrix.get_number_of_local_cols());
    std::vector<unsigned char> receive_buffer(receive_elements * sizeof(T));
    const unsigned char* receive_source
        = rank == 0 ? send_buffer.data() : nullptr;
    Parallel_Common::scatterv_bytes(receive_source,
                                    send_counts.data(),
                                    send_displacements.data(),
                                    receive_buffer.data(),
                                    checked_mpi_count(receive_buffer.size()),
                                    0,
                                    communicator);

    for (int local_row = 0; local_row < local_matrix.get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < local_matrix.get_number_of_local_cols();
             ++local_col)
        {
            const std::size_t element_index
                = static_cast<std::size_t>(local_row)
                    * static_cast<std::size_t>(local_matrix.get_number_of_local_cols())
                + static_cast<std::size_t>(local_col);
            T value;
            std::memcpy(&value, receive_buffer.data() + element_index * sizeof(T), sizeof(T));
            local_matrix.at_local(local_row, local_col) = value;
        }
    }
}
#endif
} // namespace

namespace
{
template <typename T>
void transpose_distributed(Matrix2d<T>& result, const Matrix2d<T>& source)
{
    static_cast<void>(result);
    static_cast<void>(source);
    throw std::invalid_argument(
        "distributed transpose supports only double and complex<double>");
}

template <typename T>
void matmul_distributed(Matrix2d<T>& result,
                        const Matrix2d<T>& left,
                        const Matrix2d<T>& right)
{
    static_cast<void>(result);
    static_cast<void>(left);
    static_cast<void>(right);
    throw std::invalid_argument(
        "distributed matmul supports only double and complex<double>");
}

#ifdef __MPI
void transpose_distributed(Matrix2d<double>& result,
                           const Matrix2d<double>& source)
{
    const int m = source.get_number_of_global_rows();
    const int n = source.get_number_of_global_cols();
    const double alpha = 1.0;
    const double* source_data = source.data();
    const int ia = 1;
    const int ja = 1;
    const int* source_descriptor = source.get_descriptor();
    const double beta = 0.0;
    double* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.get_descriptor();
    ScalapackConnector::tranu(n,
                              m,
                              alpha,
                              source_data,
                              ia,
                              ja,
                              source_descriptor,
                              beta,
                              result_data,
                              ic,
                              jc,
                              result_descriptor);
}

void transpose_distributed(Matrix2d<std::complex<double>>& result,
                           const Matrix2d<std::complex<double>>& source)
{
    const int m = source.get_number_of_global_rows();
    const int n = source.get_number_of_global_cols();
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double>* source_data = source.data();
    const int ia = 1;
    const int ja = 1;
    const int* source_descriptor = source.get_descriptor();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.get_descriptor();
    ScalapackConnector::tranc(n,
                              m,
                              alpha,
                              source_data,
                              ia,
                              ja,
                              source_descriptor,
                              beta,
                              result_data,
                              ic,
                              jc,
                              result_descriptor);
}

void matmul_distributed(Matrix2d<double>& result,
                        const Matrix2d<double>& left,
                        const Matrix2d<double>& right)
{
    const char transa = 'N';
    const char transb = 'N';
    const int m = left.get_number_of_global_rows();
    const int n = right.get_number_of_global_cols();
    const int k = left.get_number_of_global_cols();
    const double alpha = 1.0;
    const double* left_data = left.data();
    const int ia = 1;
    const int ja = 1;
    const int* left_descriptor = left.get_descriptor();
    const double* right_data = right.data();
    const int ib = 1;
    const int jb = 1;
    const int* right_descriptor = right.get_descriptor();
    const double beta = 0.0;
    double* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.get_descriptor();
    ScalapackConnector::gemm(transa,
                             transb,
                             m,
                             n,
                             k,
                             alpha,
                             left_data,
                             ia,
                             ja,
                             left_descriptor,
                             right_data,
                             ib,
                             jb,
                             right_descriptor,
                             beta,
                             result_data,
                             ic,
                             jc,
                             result_descriptor);
}

void matmul_distributed(Matrix2d<std::complex<double>>& result,
                        const Matrix2d<std::complex<double>>& left,
                        const Matrix2d<std::complex<double>>& right)
{
    const char transa = 'N';
    const char transb = 'N';
    const int m = left.get_number_of_global_rows();
    const int n = right.get_number_of_global_cols();
    const int k = left.get_number_of_global_cols();
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double>* left_data = left.data();
    const int ia = 1;
    const int ja = 1;
    const int* left_descriptor = left.get_descriptor();
    const std::complex<double>* right_data = right.data();
    const int ib = 1;
    const int jb = 1;
    const int* right_descriptor = right.get_descriptor();
    const std::complex<double> beta(0.0, 0.0);
    std::complex<double>* result_data = result.data();
    const int ic = 1;
    const int jc = 1;
    const int* result_descriptor = result.get_descriptor();
    ScalapackConnector::gemm(transa,
                             transb,
                             m,
                             n,
                             k,
                             alpha,
                             left_data,
                             ia,
                             ja,
                             left_descriptor,
                             right_data,
                             ib,
                             jb,
                             right_descriptor,
                             beta,
                             result_data,
                             ic,
                             jc,
                             result_descriptor);
}
#endif
} // namespace

template <typename T>
void NDArray<T, 2>::fill(const T& value)
{
    if (distribution_ == nullptr)
    {
        std::fill(storage_.begin(), storage_.end(), value);
        return;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            at_local(local_row, local_col) = value;
        }
    }
}

template <typename T>
void NDArray<T, 2>::scale(const T& value)
{
    if (distribution_ == nullptr)
    {
        for (T& element : storage_)
        {
            element *= value;
        }
        return;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            at_local(local_row, local_col) *= value;
        }
    }
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator+=(const NDArray<T, 2>& other)
{
    validate_compatible_layout(other);
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            at_local(local_row, local_col) += other.at_local(local_row, local_col);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator-=(const NDArray<T, 2>& other)
{
    validate_compatible_layout(other);
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            at_local(local_row, local_col) -= other.at_local(local_row, local_col);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::operator*=(const NDArray<T, 2>& other)
{
    validate_compatible_layout(other);
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            at_local(local_row, local_col) *= other.at_local(local_row, local_col);
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2>& NDArray<T, 2>::conjugate()
{
    if (distribution_ == nullptr)
    {
        for (int row = 0; row < rows_; ++row)
        {
            for (int col = 0; col < cols_; ++col)
            {
                at_local(row, col) = conjugate_element(at_local(row, col));
            }
        }
        return *this;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            at_local(local_row, local_col)
                = conjugate_element(at_local(local_row, local_col));
        }
    }
    return *this;
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::transpose() const
{
    if (distribution_ == nullptr || distribution_->is_serial())
    {
        // Serial degenerate path: a serial frozen distribution owns the full
        // matrix on one rank with column-major local layout; treat it like the
        // host path so distributed operators degrade to the host side.
        Matrix2d<T> result(cols_, rows_);
        for (int row = 0; row < rows_; ++row)
        {
            for (int col = 0; col < cols_; ++col)
            {
                result.at_local(col, row) = conjugate_element(at_local(row, col));
            }
        }
        return result;
    }

#ifdef __MPI
    const int result_rows = get_number_of_global_cols();
    const int result_cols = get_number_of_global_rows();
    const int block_size = distribution_->get_block_size();
    const int source_context = distribution_->get_blacs_context();
    const std::shared_ptr<const BlockCyclicDistributionInformation> result_distribution
        = make_frozen_distribution(result_rows, result_cols, block_size, source_context);
    Matrix2d<T> result(result_distribution);
    transpose_distributed(result, *this);
    return result;
#else
    throw std::invalid_argument("distributed transpose requires MPI support");
#endif
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::matmul(const NDArray<T, 2>& other) const
{
    const int left_rows = get_number_of_global_rows();
    const int shared_dimension = get_number_of_global_cols();
    const int other_shared_dimension = other.get_number_of_global_rows();
    const int result_cols = other.get_number_of_global_cols();
    if (shared_dimension != other_shared_dimension)
    {
        throw std::invalid_argument("matrix inner dimensions are incompatible");
    }

    const bool this_is_serial
        = distribution_ == nullptr || distribution_->is_serial();
    const bool other_is_serial
        = other.get_distribution() == nullptr || other.get_distribution()->is_serial();
    if (this_is_serial != other_is_serial)
    {
        throw std::invalid_argument("host and distributed matrices are incompatible");
    }

    if (this_is_serial)
    {
        // Serial degenerate path: both sides are host or serial frozen
        // distributions, so matmul degrades to the host triple loop.
        Matrix2d<T> result(left_rows, result_cols);
        for (int row = 0; row < left_rows; ++row)
        {
            for (int col = 0; col < result_cols; ++col)
            {
                T value = T();
                for (int inner = 0; inner < shared_dimension; ++inner)
                {
                    value += at_local(row, inner) * other.at_local(inner, col);
                }
                result.at_local(row, col) = value;
            }
        }
        return result;
    }

#ifdef __MPI
    const std::shared_ptr<const BlockCyclicDistributionInformation> other_distribution
        = other.get_distribution();
    if (distribution_->get_block_size() != other_distribution->get_block_size())
    {
        throw std::invalid_argument("matrix block sizes are incompatible");
    }
    if (distribution_->get_grid_rows() != other_distribution->get_grid_rows()
        || distribution_->get_grid_cols() != other_distribution->get_grid_cols())
    {
        throw std::invalid_argument("matrix process grids are incompatible");
    }
    if (distribution_->get_blacs_context() != other_distribution->get_blacs_context())
    {
        throw std::invalid_argument("matrix BLACS contexts are incompatible");
    }

    const int block_size = distribution_->get_block_size();
    const int source_context = distribution_->get_blacs_context();
    const std::shared_ptr<const BlockCyclicDistributionInformation> result_distribution
        = make_frozen_distribution(left_rows, result_cols, block_size, source_context);
    Matrix2d<T> result(result_distribution);
    matmul_distributed(result, *this, other);
    return result;
#else
    throw std::invalid_argument("distributed matmul requires MPI support");
#endif
}

template <typename T>
NDArray<double, 2> NDArray<T, 2>::real() const
{
    return extract_real_or_imaginary(*this, true);
}

template <typename T>
NDArray<double, 2> NDArray<T, 2>::imag() const
{
    return extract_real_or_imaginary(*this, false);
}

template <typename T>
void NDArray<T, 2>::real(NDArray<double, 2>& output) const
{
    validate_compatible_layout(output);
    write_real_or_imaginary(*this, output, true);
}

template <typename T>
void NDArray<T, 2>::imag(NDArray<double, 2>& output) const
{
    validate_compatible_layout(output);
    write_real_or_imaginary(*this, output, false);
}

template <typename T>
T NDArray<T, 2>::sum() const
{
    T value = T();
    if (distribution_ == nullptr)
    {
        for (const T& element : storage_)
        {
            value += element;
        }
        return value;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            value += at_local(local_row, local_col);
        }
    }
#ifdef __MPI
    if (!distribution_->is_serial())
    {
        const MPI_Comm communicator = distribution_->get_mpi_communicator();
        if (communicator != MPI_COMM_NULL)
        {
            Parallel_Common::allreduce_sum(value, communicator);
        }
    }
#endif
    return value;
}

template <typename T>
double NDArray<T, 2>::sum_squares() const
{
    double value = 0.0;
    if (distribution_ == nullptr)
    {
        for (int row = 0; row < rows_; ++row)
        {
            for (int col = 0; col < cols_; ++col)
            {
                value += squared_element(at_local(row, col));
            }
        }
        return value;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            value += squared_element(at_local(local_row, local_col));
        }
    }
#ifdef __MPI
    if (!distribution_->is_serial())
    {
        const MPI_Comm communicator = distribution_->get_mpi_communicator();
        if (communicator != MPI_COMM_NULL)
        {
            Parallel_Common::allreduce_sum(value, communicator);
        }
    }
#endif
    return value;
}

template <typename T>
T NDArray<T, 2>::dot(const NDArray<T, 2>& other) const
{
    validate_compatible_layout(other);
    T value = T();
    if (distribution_ == nullptr)
    {
        for (int row = 0; row < rows_; ++row)
        {
            for (int col = 0; col < cols_; ++col)
            {
                value += conjugate_element(at_local(row, col))
                    * other.at_local(row, col);
            }
        }
        return value;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            value += conjugate_element(at_local(local_row, local_col))
                * other.at_local(local_row, local_col);
        }
    }
#ifdef __MPI
    if (!distribution_->is_serial())
    {
        const MPI_Comm communicator = distribution_->get_mpi_communicator();
        if (communicator != MPI_COMM_NULL)
        {
            Parallel_Common::allreduce_sum(value, communicator);
        }
    }
#endif
    return value;
}

template <typename T>
double NDArray<T, 2>::norm2() const
{
    return std::sqrt(sum_squares());
}

template <typename T>
double NDArray<T, 2>::max_abs() const
{
    double value = 0.0;
    if (distribution_ == nullptr)
    {
        for (int row = 0; row < rows_; ++row)
        {
            for (int col = 0; col < cols_; ++col)
            {
                value = std::max(value,
                                 static_cast<double>(std::abs(at_local(row, col))));
            }
        }
        return value;
    }
    for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
    {
        for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
        {
            value = std::max(
                value,
                static_cast<double>(std::abs(at_local(local_row, local_col))));
        }
    }
#ifdef __MPI
    if (!distribution_->is_serial())
    {
        const MPI_Comm communicator = distribution_->get_mpi_communicator();
        if (communicator != MPI_COMM_NULL)
        {
            Parallel_Common::allreduce_max(value, communicator);
        }
    }
#endif
    return value;
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::gather() const
{
    if (distribution_ == nullptr)
    {
        return *this;
    }
    if (distribution_->is_serial())
    {
        Matrix2d<T> result(get_number_of_global_rows(), get_number_of_global_cols());
        for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
        {
            const int global_row = distribution_->local_to_global_row(local_row);
            for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
            {
                const int global_col = distribution_->local_to_global_col(local_col);
                result.at_local(global_row, global_col) = at_local(local_row, local_col);
            }
        }
        return result;
    }
#ifdef __MPI
    const MPI_Comm communicator = distribution_->get_mpi_communicator();
    if (communicator == MPI_COMM_NULL)
    {
        throw std::invalid_argument("gather requires a valid MPI communicator");
    }
    const std::vector<std::vector<unsigned char>> payloads
        = allgather_payloads(*this, communicator);
    return rebuild_global_matrix<T>(payloads);
#else
    return *this;
#endif
}

template <typename T>
NDArray<T, 2> NDArray<T, 2>::gather_on_root(const int root_rank) const
{
    if (distribution_ == nullptr)
    {
        if (root_rank == 0)
        {
            return *this;
        }
        return NDArray(EmptyMatrixTag());
    }
    if (distribution_->is_serial())
    {
        if (root_rank != 0)
        {
            return NDArray(EmptyMatrixTag());
        }
        Matrix2d<T> result(get_number_of_global_rows(), get_number_of_global_cols());
        for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
        {
            const int global_row = distribution_->local_to_global_row(local_row);
            for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
            {
                const int global_col = distribution_->local_to_global_col(local_col);
                result.at_local(global_row, global_col) = at_local(local_row, local_col);
            }
        }
        return result;
    }
#ifdef __MPI
    const MPI_Comm communicator = distribution_->get_mpi_communicator();
    if (communicator == MPI_COMM_NULL)
    {
        throw std::invalid_argument("gather_on_root requires a valid MPI communicator");
    }
    std::unique_ptr<Matrix2d<T>> gathered
        = gather_payloads_on_root(*this, root_rank, communicator);
    if (gathered)
    {
        return *gathered;
    }
    return NDArray(EmptyMatrixTag());
#else
    static_cast<void>(root_rank);
    return NDArray(EmptyMatrixTag());
#endif
}

template <typename T>
void NDArray<T, 2>::scatter(const NDArray<T, 2>& global_matrix)
{
    if (distribution_ == nullptr)
    {
        if (global_matrix.is_distributed())
        {
            throw std::invalid_argument("host scatter requires a host global matrix");
        }
        if (rows_ != global_matrix.get_number_of_global_rows()
            || cols_ != global_matrix.get_number_of_global_cols())
        {
            throw std::invalid_argument("scatter global matrix shape does not match matrix");
        }
        *this = global_matrix;
        return;
    }

    if (distribution_->is_serial())
    {
        if (global_matrix.is_distributed())
        {
            throw std::invalid_argument("scatter requires a host global matrix");
        }
        if (get_number_of_global_rows() != global_matrix.get_number_of_global_rows()
            || get_number_of_global_cols() != global_matrix.get_number_of_global_cols())
        {
            throw std::invalid_argument("scatter global matrix shape does not match distribution");
        }
        for (int local_row = 0; local_row < get_number_of_local_rows(); ++local_row)
        {
            for (int local_col = 0; local_col < get_number_of_local_cols(); ++local_col)
            {
                at_local(local_row, local_col)
                    = global_matrix.at_local(distribution_->local_to_global_row(local_row),
                                             distribution_->local_to_global_col(local_col));
            }
        }
        return;
    }

#ifdef __MPI
    scatter_from_host(global_matrix, *this);
#else
    throw std::invalid_argument("active distributed scatter requires MPI support");
#endif
}

template class NDArray<int, 2>;
template class NDArray<float, 2>;
template class NDArray<double, 2>;
template class NDArray<std::complex<float>, 2>;
template class NDArray<std::complex<double>, 2>;
