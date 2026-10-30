#ifdef NDEBUG
#undef NDEBUG
#endif

#include "source_base/module_container/NDArray/distributed_ndarray.h"
#include "source_base/parallel_common.h"

#include "gtest/gtest.h"

#include <cstddef>
#include <memory>
#include <type_traits>

static_assert(!std::is_copy_constructible<DistributedNDArray<double>>::value,
              "DistributedNDArray must not be copy constructible");
static_assert(!std::is_copy_assignable<DistributedNDArray<double>>::value,
              "DistributedNDArray must not be copy assignable");
static_assert(std::is_nothrow_move_constructible<DistributedNDArray<double>>::value,
              "DistributedNDArray move construction must be noexcept");
static_assert(std::is_nothrow_move_assignable<DistributedNDArray<double>>::value,
              "DistributedNDArray move assignment must be noexcept");

namespace
{
#ifdef __MPI
constexpr int kDescriptorMatrixTypeIndex = 0;
constexpr int kDescriptorBlacsContextIndex = 1;
constexpr int kDescriptorGlobalRowsIndex = 2;
constexpr int kDescriptorGlobalColsIndex = 3;
constexpr int kDescriptorBlockRowsIndex = 4;
constexpr int kDescriptorBlockColsIndex = 5;
constexpr int kDescriptorProcessRowsIndex = 6;
constexpr int kDescriptorProcessColsIndex = 7;
constexpr int kDescriptorLeadingDimIndex = 8;

void choose_process_grid(int& number_of_row_processors, int& number_of_column_processors)
{
    const int process_count = Parallel_Common::comm_size(MPI_COMM_WORLD);
    number_of_row_processors = 1;
    number_of_column_processors = process_count;
    for (int candidate = 1; candidate * candidate <= process_count; ++candidate)
    {
        if (process_count % candidate == 0)
        {
            number_of_row_processors = candidate;
            number_of_column_processors = process_count / candidate;
        }
    }
}
#endif

std::shared_ptr<BlockCyclicDistributionInformation> make_serial_distribution(
    const int global_rows,
    const int global_cols)
{
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->initialize_serial(global_rows, global_cols);

#ifdef __MPI
    int* const descriptor = distribution->get_mutable_descriptor();
    descriptor[kDescriptorMatrixTypeIndex] = 1;
    descriptor[kDescriptorBlacsContextIndex] = -1;
    descriptor[kDescriptorGlobalRowsIndex] = global_rows;
    descriptor[kDescriptorGlobalColsIndex] = global_cols;
    descriptor[kDescriptorBlockRowsIndex] = 1;
    descriptor[kDescriptorBlockColsIndex] = 1;
    descriptor[kDescriptorProcessRowsIndex] = 0;
    descriptor[kDescriptorProcessColsIndex] = 0;
    descriptor[kDescriptorLeadingDimIndex] = global_rows > 1 ? global_rows : 1;
#endif

    return distribution;
}

void check_linearization(const int global_rows,
                         const int global_cols,
                         const StorageMajorOrder major_order)
{
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    const int expected_leading_dimension = global_rows > 1 ? global_rows : 1;

    DistributedNDArray<double> array(distribution, major_order, expected_leading_dimension);

    EXPECT_EQ(array.get_distribution().get(), distribution.get());
    EXPECT_EQ(array.get_major_order(), major_order);
    EXPECT_EQ(array.get_leading_dimension(), expected_leading_dimension);
    EXPECT_EQ(array.get_number_of_local_rows(), global_rows);
    EXPECT_EQ(array.get_number_of_local_cols(), global_cols);
    EXPECT_EQ(array.get_local_block().data(), array.data());

    std::size_t leading_dimension = static_cast<std::size_t>(array.get_leading_dimension());
#ifdef __MPI
    const int* const descriptor = distribution->get_descriptor();
    EXPECT_EQ(array.get_leading_dimension(), descriptor[kDescriptorLeadingDimIndex]);
    leading_dimension = static_cast<std::size_t>(descriptor[kDescriptorLeadingDimIndex]);
#endif

    for (int local_row = 0; local_row < global_rows; ++local_row)
    {
        for (int local_col = 0; local_col < global_cols; ++local_col)
        {
            std::size_t expected_flat_index = 0;
            if (major_order == StorageMajorOrder::ColumnMajor)
            {
                expected_flat_index = static_cast<std::size_t>(local_col) * leading_dimension
                    + static_cast<std::size_t>(local_row);
            }
            else
            {
                expected_flat_index = static_cast<std::size_t>(local_row) * leading_dimension
                    + static_cast<std::size_t>(local_col);
            }

            const double* const expected_address = array.data() + expected_flat_index;
            EXPECT_EQ(&array.at_local(local_row, local_col), expected_address);

            const double element_value = static_cast<double>(local_row * global_cols + local_col);
            array.at_local(local_row, local_col) = element_value;
            EXPECT_EQ(array.at_local(local_row, local_col), element_value);

            const int mapped_local_row = distribution->global_to_local_row(local_row);
            const int mapped_local_col = distribution->global_to_local_col(local_col);
            ASSERT_NE(mapped_local_row, -1);
            ASSERT_NE(mapped_local_col, -1);
            EXPECT_EQ(&array.at_global(local_row, local_col),
                      &array.at_local(mapped_local_row, mapped_local_col));
            EXPECT_EQ(array.at_global(local_row, local_col),
                      array.at_local(mapped_local_row, mapped_local_col));

            const bool expected_owner = distribution->owns_index(local_row, local_col);
            EXPECT_EQ(array.owns_global_index(local_row, local_col), expected_owner);
        }
    }
}

#ifdef __MPI
void check_mpi_linearization(const int global_rows,
                             const int global_cols,
                             const int block_size,
                             const StorageMajorOrder major_order)
{
    int number_of_row_processors = 0;
    int number_of_column_processors = 0;
    choose_process_grid(number_of_row_processors, number_of_column_processors);

    const std::shared_ptr<const BlockCyclicDistributionInformation> distribution
        = make_frozen_distribution(global_rows,
                                   global_cols,
                                   number_of_row_processors,
                                   number_of_column_processors,
                                   block_size,
                                   MPI_COMM_WORLD);
    ASSERT_NE(distribution, nullptr);

    const int leading_dimension = distribution->get_leading_dimension();
    DistributedNDArray<double> array(distribution, major_order, leading_dimension);
    EXPECT_EQ(array.get_leading_dimension(), leading_dimension);
    EXPECT_EQ(array.get_distribution().get(), distribution.get());
    EXPECT_EQ(array.get_number_of_local_rows(), distribution->get_number_of_local_rows());
    EXPECT_EQ(array.get_number_of_local_cols(), distribution->get_number_of_local_cols());

    for (int global_row = 0; global_row < global_rows; ++global_row)
    {
        for (int global_col = 0; global_col < global_cols; ++global_col)
        {
            const bool owns_index = distribution->owns_index(global_row, global_col);
            const bool array_owns_index = array.owns_global_index(global_row, global_col);
            EXPECT_EQ(array_owns_index, owns_index);

            const int local_row = distribution->global_to_local_row(global_row);
            const int local_col = distribution->global_to_local_col(global_col);
            const bool local_owns_row = local_row != -1;
            const bool local_owns_col = local_col != -1;
            const int row_owner = (global_row / block_size) % number_of_row_processors;
            const int column_owner = (global_col / block_size) % number_of_column_processors;
            const bool expected_local_row = row_owner == distribution->get_coord_row();
            const bool expected_local_col = column_owner == distribution->get_coord_col();
            EXPECT_EQ(local_owns_row, expected_local_row);
            EXPECT_EQ(local_owns_col, expected_local_col);

            if (owns_index)
            {
                ASSERT_NE(local_row, -1);
                ASSERT_NE(local_col, -1);
                double& global_element = array.at_global(global_row, global_col);
                const double& local_element = array.at_local(local_row, local_col);
                EXPECT_EQ(&global_element, &local_element);

                const double element_value
                    = static_cast<double>(global_row * global_cols + global_col);
                global_element = element_value;
                EXPECT_EQ(array.at_local(local_row, local_col), element_value);
            }
        }
    }
}
#endif
} // namespace

TEST(DistributedNDArrayTest, ColumnMajorLinearization)
{
    const int global_rows = 8;
    const int global_cols = 6;
    check_linearization(global_rows, global_cols, StorageMajorOrder::ColumnMajor);
}

TEST(DistributedNDArrayTest, RowMajorLinearization)
{
    const int global_rows = 5;
    const int global_cols = 3;
    check_linearization(global_rows, global_cols, StorageMajorOrder::RowMajor);
}

TEST(DistributedNDArrayTest, AtGlobalAssertsOutsideGlobalBounds)
{
    const int global_rows = 4;
    const int global_cols = 3;
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    const int leading_dimension = global_rows > 1 ? global_rows : 1;
    DistributedNDArray<double> array(distribution, StorageMajorOrder::ColumnMajor, leading_dimension);

    EXPECT_DEATH(array.at_global(global_rows, 0), "");
}

#ifdef __MPI
TEST(DistributedNDArrayTest, MpiAtLocalAtGlobal8x6)
{
    check_mpi_linearization(8, 6, 2, StorageMajorOrder::ColumnMajor);
}

TEST(DistributedNDArrayTest, MpiAtLocalAtGlobal5x3)
{
    check_mpi_linearization(5, 3, 2, StorageMajorOrder::RowMajor);
}

TEST(DistributedNDArrayTest, RejectsMismatchedDescriptorLeadingDimension)
{
    const int global_rows = 4;
    const int global_cols = 3;
    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = make_serial_distribution(global_rows, global_cols);
    int* const descriptor = distribution->get_mutable_descriptor();
    descriptor[kDescriptorLeadingDimIndex] += 1;
    const int expected_leading_dimension = global_rows > 1 ? global_rows : 1;

    EXPECT_DEATH(
        {
            DistributedNDArray<double> invalid_array(
                distribution,
                StorageMajorOrder::ColumnMajor,
                expected_leading_dimension);
            static_cast<void>(invalid_array);
        },
        "");
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
