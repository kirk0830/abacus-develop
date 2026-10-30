#ifndef ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_H
#define ABACUS_SOURCE_SOURCE_BASE_MATRIX2D_H

#include "source_base/block_cyclic_distribution_information.h"
#include "source_base/ndarray.h"

#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

template <typename T>
class NDArray<T, 2> final
{
  public:
    NDArray() = delete;

    NDArray(const int rows, const int cols)
        : rows_(rows),
          cols_(cols)
    {
        if (rows_ <= 0)
        {
            throw std::invalid_argument("rows must be positive");
        }
        if (cols_ <= 0)
        {
            throw std::invalid_argument("cols must be positive");
        }
        const std::size_t number_of_elements
            = static_cast<std::size_t>(rows_) * static_cast<std::size_t>(cols_);
        storage_.resize(number_of_elements);
    }

#ifdef __MPI
    NDArray(const int global_rows,
            const int global_cols,
            const MPI_Comm communicator)
    {
        const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
            = make_frozen_distribution(global_rows, global_cols, 1, communicator);
        initialize_distributed_storage(distribution);
    }

    NDArray(const int global_rows,
            const int global_cols,
            const int block_size,
            const MPI_Comm communicator)
    {
        const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
            = make_frozen_distribution(global_rows, global_cols, block_size, communicator);
        initialize_distributed_storage(distribution);
    }

    NDArray(const int global_rows,
            const int global_cols,
            const int number_of_row_processors,
            const int number_of_column_processors,
            const int block_size,
            const MPI_Comm communicator)
    {
        const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
            = make_frozen_distribution(global_rows,
                                       global_cols,
                                       number_of_row_processors,
                                       number_of_column_processors,
                                       block_size,
                                       communicator);
        initialize_distributed_storage(distribution);
    }
#endif

    explicit NDArray(std::shared_ptr<const BlockCyclicDistributionInformation> distribution)
        : distribution_(std::move(distribution))
    {
        if (distribution_ == nullptr)
        {
            throw std::invalid_argument("distribution must not be null");
        }
        initialize_distributed_storage();
    }

    NDArray(const NDArray&) = default;
    NDArray(NDArray&&) noexcept = default;
    ~NDArray() = default;

    NDArray& operator=(const NDArray&) = default;
    NDArray& operator=(NDArray&&) noexcept = default;

    void fill(const T& value);
    void scale(const T& value);
    NDArray& operator+=(const NDArray& other);
    NDArray& operator-=(const NDArray& other);
    NDArray& operator*=(const NDArray& other);

    NDArray& conjugate();

    NDArray<double, 2> real() const;
    NDArray<double, 2> imag() const;
    void real(NDArray<double, 2>& output) const;
    void imag(NDArray<double, 2>& output) const;

    T sum() const;
    double sum_squares() const;
    T dot(const NDArray& other) const;
    double norm2() const;
    double max_abs() const;

    NDArray gather() const;
    NDArray gather_on_root(const int root_rank) const;
    void scatter(const NDArray& global_matrix);

    T& at_global(const int global_row, const int global_col)
    {
        assert(global_row >= 0);
        assert(global_row < get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < get_number_of_global_cols());

        if (distribution_ == nullptr)
        {
            return at_local(global_row, global_col);
        }

        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        const bool owns_index = distribution_->owns_index(global_row_index, global_col_index);
        assert(owns_index);

        const int local_row = distribution_->global_to_local_row(global_row);
        const int local_col = distribution_->global_to_local_col(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    const T& at_global(const int global_row, const int global_col) const
    {
        assert(global_row >= 0);
        assert(global_row < get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < get_number_of_global_cols());

        if (distribution_ == nullptr)
        {
            return at_local(global_row, global_col);
        }

        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        const bool owns_index = distribution_->owns_index(global_row_index, global_col_index);
        assert(owns_index);

        const int local_row = distribution_->global_to_local_row(global_row);
        const int local_col = distribution_->global_to_local_col(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    T& at_local(const int local_row, const int local_col)
    {
        assert(local_row >= 0);
        assert(local_row < get_number_of_local_rows());
        assert(local_col >= 0);
        assert(local_col < get_number_of_local_cols());

        std::size_t storage_index = 0;
        if (distribution_ == nullptr)
        {
            storage_index = static_cast<std::size_t>(local_row) * static_cast<std::size_t>(cols_)
                + static_cast<std::size_t>(local_col);
        }
        else
        {
            const int leading_dimension = distribution_->get_leading_dimension();
            storage_index = static_cast<std::size_t>(local_col) * static_cast<std::size_t>(leading_dimension)
                + static_cast<std::size_t>(local_row);
        }
        assert(storage_index < storage_.size());
        return storage_[storage_index];
    }

    const T& at_local(const int local_row, const int local_col) const
    {
        assert(local_row >= 0);
        assert(local_row < get_number_of_local_rows());
        assert(local_col >= 0);
        assert(local_col < get_number_of_local_cols());

        std::size_t storage_index = 0;
        if (distribution_ == nullptr)
        {
            storage_index = static_cast<std::size_t>(local_row) * static_cast<std::size_t>(cols_)
                + static_cast<std::size_t>(local_col);
        }
        else
        {
            const int leading_dimension = distribution_->get_leading_dimension();
            storage_index = static_cast<std::size_t>(local_col) * static_cast<std::size_t>(leading_dimension)
                + static_cast<std::size_t>(local_row);
        }
        assert(storage_index < storage_.size());
        return storage_[storage_index];
    }

    bool owns_global_index(const int global_row, const int global_col) const
    {
        assert(global_row >= 0);
        assert(global_row < get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < get_number_of_global_cols());
        if (distribution_ == nullptr)
        {
            return true;
        }
        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        return distribution_->owns_index(global_row_index, global_col_index);
    }

    T* data()
    {
        return storage_.data();
    }

    const T* data() const
    {
        return storage_.data();
    }

    int get_number_of_global_rows() const
    {
        if (distribution_ == nullptr)
        {
            return rows_;
        }
        return distribution_->get_number_of_global_rows();
    }

    int get_number_of_global_cols() const
    {
        if (distribution_ == nullptr)
        {
            return cols_;
        }
        return distribution_->get_number_of_global_cols();
    }

    int get_number_of_local_rows() const
    {
        if (distribution_ == nullptr)
        {
            return rows_;
        }
        return distribution_->get_number_of_local_rows();
    }

    int get_number_of_local_cols() const
    {
        if (distribution_ == nullptr)
        {
            return cols_;
        }
        return distribution_->get_number_of_local_cols();
    }

    int64_t get_number_of_local_elements() const
    {
        if (distribution_ != nullptr)
        {
            return distribution_->get_number_of_local_elements();
        }
        const int64_t number_of_rows = static_cast<int64_t>(rows_);
        const int64_t number_of_cols = static_cast<int64_t>(cols_);
        return number_of_rows * number_of_cols;
    }

    int get_grid_rows() const
    {
        if (distribution_ == nullptr)
        {
            return 1;
        }
        return distribution_->get_grid_rows();
    }

    int get_grid_cols() const
    {
        if (distribution_ == nullptr)
        {
            return 1;
        }
        return distribution_->get_grid_cols();
    }

    int get_coord_row() const
    {
        if (distribution_ == nullptr)
        {
            return 0;
        }
        return distribution_->get_coord_row();
    }

    int get_coord_col() const
    {
        if (distribution_ == nullptr)
        {
            return 0;
        }
        return distribution_->get_coord_col();
    }

    int get_leading_dimension() const
    {
        if (distribution_ == nullptr)
        {
            return cols_;
        }
        return distribution_->get_leading_dimension();
    }

    const int* get_descriptor() const
    {
#ifdef __MPI
        if (distribution_ != nullptr)
        {
            return distribution_->get_descriptor();
        }
#endif
        return nullptr;
    }

    int get_blacs_context() const
    {
#ifdef __MPI
        if (distribution_ != nullptr)
        {
            return distribution_->get_blacs_context();
        }
#endif
        return -1;
    }

    const std::shared_ptr<const BlockCyclicDistributionInformation>& get_distribution() const
    {
        return distribution_;
    }

    bool is_distributed() const
    {
        return distribution_ != nullptr;
    }

  private:
    struct EmptyMatrixTag
    {
    };

    explicit NDArray(EmptyMatrixTag)
    {
    }

    template <typename U>
    void validate_compatible_layout(const NDArray<U, 2>& other) const
    {
        if (get_number_of_global_rows() != other.get_number_of_global_rows()
            || get_number_of_global_cols() != other.get_number_of_global_cols())
        {
            throw std::invalid_argument("matrix global shapes are incompatible");
        }

        const bool this_is_distributed = distribution_ != nullptr;
        const bool other_is_distributed = other.get_distribution() != nullptr;
        if (this_is_distributed != other_is_distributed)
        {
            throw std::invalid_argument("host and distributed matrices are incompatible");
        }

        if (!this_is_distributed)
        {
            return;
        }

        if (distribution_->get_block_size() != other.get_distribution()->get_block_size())
        {
            throw std::invalid_argument("matrix block sizes are incompatible");
        }
        if (distribution_->get_grid_rows() != other.get_distribution()->get_grid_rows()
            || distribution_->get_grid_cols() != other.get_distribution()->get_grid_cols())
        {
            throw std::invalid_argument("matrix process grids are incompatible");
        }
        if (distribution_ != other.get_distribution())
        {
            throw std::invalid_argument("matrix distribution handles are incompatible");
        }
    }

    void initialize_distributed_storage(
        const std::shared_ptr<const BlockCyclicDistributionInformation>& distribution)
    {
        if (distribution == nullptr)
        {
            throw std::invalid_argument("distribution must not be null");
        }
        distribution_ = distribution;
        initialize_distributed_storage();
    }

    void initialize_distributed_storage()
    {
        rows_ = distribution_->get_number_of_local_rows();
        cols_ = distribution_->get_number_of_local_cols();
        assert(rows_ >= 0);
        assert(cols_ >= 0);

        const int leading_dimension = std::max(rows_, 1);
        const int row_major_leading_dimension = std::max(cols_, leading_dimension);
        const std::size_t column_major_storage_elements
            = static_cast<std::size_t>(leading_dimension) * static_cast<std::size_t>(cols_);
        const std::size_t row_major_storage_elements
            = static_cast<std::size_t>(rows_) * static_cast<std::size_t>(row_major_leading_dimension);
        const std::size_t number_of_storage_elements
            = std::max(column_major_storage_elements, row_major_storage_elements);
        storage_.resize(number_of_storage_elements);
    }

    int rows_ = 0;
    int cols_ = 0;
    std::vector<T> storage_;
    std::shared_ptr<const BlockCyclicDistributionInformation> distribution_;
};

template <typename T>
using Matrix2d = NDArray<T, 2>;

static_assert(std::is_same<Matrix2d<int>, NDArray<int, 2>>::value,
              "Matrix2d<T> and NDArray<T, 2> must be the same type.");
static_assert(!std::is_same<NDArray<int>, NDArray<int, 2>>::value,
              "NDArray<T> and NDArray<T, 2> must be distinct types.");
static_assert(!std::is_base_of<NDArray<int>, NDArray<int, 2>>::value,
              "The rank-two specialization must compose storage rather than inherit.");

#endif
