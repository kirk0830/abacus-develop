#ifndef ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DISTRIBUTED_NDARRAY_H
#define ABACUS_SOURCE_SOURCE_BASE_MODULE_CONTAINER_NDARRAY_DISTRIBUTED_NDARRAY_H

#include "source_base/matrix2d.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>

enum class StorageMajorOrder
{
    RowMajor,
    ColumnMajor
};

template <typename T, typename Tdata = double>
class DistributedNDArray
{
  public:
    DistributedNDArray(std::shared_ptr<const BlockCyclicDistributionInformation> distribution,
                       const StorageMajorOrder major_order,
                       const int leading_dimension)
        : local_block_(std::move(distribution)),
          major_order_(major_order)
    {
        const int matrix_leading_dimension = local_block_.get_leading_dimension();
        assert(leading_dimension == matrix_leading_dimension);
#ifdef __MPI
        const int* const descriptor = local_block_.get_descriptor();
        assert(descriptor != nullptr);
        constexpr int kDescriptorLeadingDimIndex = 8;
        assert(leading_dimension == descriptor[kDescriptorLeadingDimIndex]);
#endif
    }

    ~DistributedNDArray() = default;

    DistributedNDArray(const DistributedNDArray&) = delete;
    DistributedNDArray& operator=(const DistributedNDArray&) = delete;

    DistributedNDArray(DistributedNDArray&&) noexcept = default;
    DistributedNDArray& operator=(DistributedNDArray&&) noexcept = default;

    T& at_local(const int local_row, const int local_col)
    {
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            return local_block_.at_local(local_row, local_col);
        }
        return at_local_row_major(local_row, local_col);
    }

    const T& at_local(const int local_row, const int local_col) const
    {
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            return local_block_.at_local(local_row, local_col);
        }
        return at_local_row_major(local_row, local_col);
    }

    T& at_global(const int global_row, const int global_col)
    {
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            return local_block_.at_global(global_row, global_col);
        }

        assert(global_row >= 0);
        assert(global_row < local_block_.get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < local_block_.get_number_of_global_cols());
        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        const std::shared_ptr<const BlockCyclicDistributionInformation>& distribution
            = local_block_.get_distribution();
        assert(distribution->owns_index(global_row_index, global_col_index));

        const int local_row = distribution->global_to_local_row(global_row);
        const int local_col = distribution->global_to_local_col(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    const T& at_global(const int global_row, const int global_col) const
    {
        if (major_order_ == StorageMajorOrder::ColumnMajor)
        {
            return local_block_.at_global(global_row, global_col);
        }

        assert(global_row >= 0);
        assert(global_row < local_block_.get_number_of_global_rows());
        assert(global_col >= 0);
        assert(global_col < local_block_.get_number_of_global_cols());
        const std::size_t global_row_index = static_cast<std::size_t>(global_row);
        const std::size_t global_col_index = static_cast<std::size_t>(global_col);
        const std::shared_ptr<const BlockCyclicDistributionInformation>& distribution
            = local_block_.get_distribution();
        assert(distribution->owns_index(global_row_index, global_col_index));

        const int local_row = distribution->global_to_local_row(global_row);
        const int local_col = distribution->global_to_local_col(global_col);
        assert(local_row != -1);
        assert(local_col != -1);
        return at_local(local_row, local_col);
    }

    bool owns_global_index(const int global_row, const int global_col) const
    {
        return local_block_.owns_global_index(global_row, global_col);
    }

    T* data()
    {
        return local_block_.data();
    }

    const T* data() const
    {
        return local_block_.data();
    }

    const std::shared_ptr<const BlockCyclicDistributionInformation>& get_distribution() const
    {
        return local_block_.get_distribution();
    }

    Matrix2d<T>& get_local_block()
    {
        return local_block_;
    }

    const Matrix2d<T>& get_local_block() const
    {
        return local_block_;
    }

    StorageMajorOrder get_major_order() const
    {
        return major_order_;
    }

    int get_leading_dimension() const
    {
        return local_block_.get_leading_dimension();
    }

    int get_number_of_local_rows() const
    {
        return local_block_.get_number_of_local_rows();
    }

    int get_number_of_local_cols() const
    {
        return local_block_.get_number_of_local_cols();
    }

  private:
    T& at_local_row_major(const int local_row, const int local_col)
    {
        const std::size_t storage_index = row_major_index(local_row, local_col);
        return data()[storage_index];
    }

    const T& at_local_row_major(const int local_row, const int local_col) const
    {
        const std::size_t storage_index = row_major_index(local_row, local_col);
        return data()[storage_index];
    }

    std::size_t row_major_index(const int local_row, const int local_col) const
    {
        assert(local_row >= 0);
        assert(local_row < local_block_.get_number_of_local_rows());
        assert(local_col >= 0);
        assert(local_col < local_block_.get_number_of_local_cols());
        const std::size_t row = static_cast<std::size_t>(local_row);
        const std::size_t leading_dimension
            = static_cast<std::size_t>(local_block_.get_leading_dimension());
        const std::size_t column = static_cast<std::size_t>(local_col);
        return row * leading_dimension + column;
    }

    Matrix2d<T> local_block_;
    StorageMajorOrder major_order_;
};

#endif
