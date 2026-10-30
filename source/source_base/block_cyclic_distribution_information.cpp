#include "source_base/block_cyclic_distribution_information.h"

#include "source_base/module_external/blacs_connector.h"
#include "source_base/module_external/scalapack_connector.h"
#include "source_base/parallel_common.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace
{
void validate_global_shape(const int global_rows, const int global_cols)
{
    if (global_rows <= 0)
    {
        throw std::invalid_argument("global_rows must be positive");
    }
    if (global_cols <= 0)
    {
        throw std::invalid_argument("global_cols must be positive");
    }
}

#ifdef __MPI
void validate_block_size(const int block_size)
{
    if (block_size <= 0)
    {
        throw std::invalid_argument("block_size must be positive");
    }
}
#endif
} // namespace

BlockCyclicDistributionInformation::~BlockCyclicDistributionInformation()
{
#ifdef __MPI
    release_blacs_grid();
#endif
}

BlockCyclicDistributionInformation::BlockCyclicDistributionInformation(
    BlockCyclicDistributionInformation&& rhs) noexcept
{
    *this = std::move(rhs);
}

BlockCyclicDistributionInformation& BlockCyclicDistributionInformation::operator=(
    BlockCyclicDistributionInformation&& rhs) noexcept
{
    if (this != &rhs)
    {
#ifdef __MPI
        release_blacs_grid();
        owns_blacs_context_ = rhs.owns_blacs_context_;
        rhs.owns_blacs_context_ = false;
        blacs_context_ = rhs.blacs_context_;
        rhs.blacs_context_ = kUnownedBlacsContext;
#endif
        number_of_local_rows_ = rhs.number_of_local_rows_;
        number_of_local_cols_ = rhs.number_of_local_cols_;
        number_of_local_elements_ = rhs.number_of_local_elements_;
        block_size_ = rhs.block_size_;
        grid_rows_ = rhs.grid_rows_;
        grid_cols_ = rhs.grid_cols_;
        coord_[0] = rhs.coord_[0];
        coord_[1] = rhs.coord_[1];
        is_serial_ = rhs.is_serial_;
        global_to_local_row_ = std::move(rhs.global_to_local_row_);
        global_to_local_col_ = std::move(rhs.global_to_local_col_);
        local_to_global_row_ = std::move(rhs.local_to_global_row_);
        local_to_global_col_ = std::move(rhs.local_to_global_col_);
#ifdef __MPI
        std::copy(rhs.descriptor_.begin(), rhs.descriptor_.end(), descriptor_.begin());
#endif
    }
    return *this;
}

bool BlockCyclicDistributionInformation::owns_index(const std::size_t iw1_all,
                                                    const std::size_t iw2_all) const
{
    return global2local_row(iw1_all) != -1 && global2local_col(iw2_all) != -1;
}

int BlockCyclicDistributionInformation::owner_of_index(const std::size_t iw1_all,
                                                       const std::size_t iw2_all) const
{
    assert(iw1_all < get_global_row_size() && iw2_all < get_global_col_size());
    if (is_serial_)
    {
        return 0;
    }
    int proc_row = (iw1_all / block_size_) % grid_rows_;
    int proc_col = (iw2_all / block_size_) % grid_cols_;
    return proc_row * grid_cols_ + proc_col;
}

int BlockCyclicDistributionInformation::get_number_of_global_rows() const
{
    if (!is_serial_)
    {
#ifdef __MPI
        return descriptor_[kDescriptorGlobalRowsIndex];
#endif
    }
    return number_of_local_rows_;
}

int BlockCyclicDistributionInformation::get_number_of_global_cols() const
{
    if (!is_serial_)
    {
#ifdef __MPI
        return descriptor_[kDescriptorGlobalColsIndex];
#endif
    }
    return number_of_local_cols_;
}

int BlockCyclicDistributionInformation::get_leading_dimension() const
{
#ifdef __MPI
    if (!is_serial_ && blacs_context_ >= 0)
    {
        return descriptor_[kDescriptorLeadingDimIndex];
    }
#endif
    return std::max(number_of_local_rows_, 1);
}

#ifdef __MPI
MPI_Comm BlockCyclicDistributionInformation::get_mpi_communicator() const
{
    // it is an error to call blacs_get with an invalid BLACS context
    if (blacs_context_ < 0)
    {
        return MPI_COMM_NULL;
    }

    int sys_ctxt = 0;
    Cblacs_get(blacs_context_, kBlacsWhatSystemContext, &sys_ctxt);
    // blacs_get with "what" = 10 takes a BLACS context and returns the index
    // of the associated system context (MPI communicator) that can be used by
    // blacs2sys_handle to get the MPI communicator.
    return Cblacs2sys_handle(sys_ctxt);
}

void BlockCyclicDistributionInformation::initialize_blacs_grid(const MPI_Comm comm,
                                                               const bool swap_grid)
{
    // determine the number of rows and columns of the process grid
    // by factorizing n = p * q such that p, q are closest and p <= q
    const int num_proc = Parallel_Common::comm_size(comm);
    grid_rows_ = static_cast<int>(std::sqrt(num_proc + 0.5));
    while (grid_cols_ = num_proc / grid_rows_, grid_rows_ * grid_cols_ != num_proc)
    {
        --grid_rows_;
    }

    if (swap_grid)
    {
        std::swap(grid_rows_, grid_cols_);
    }

    // initialize the BLACS grid accordingly
    blacs_context_ = Csys2blacs_handle(comm);
    char order = 'R'; // row-major
    Cblacs_gridinit(&blacs_context_, &order, grid_rows_, grid_cols_);
    Cblacs_gridinfo(blacs_context_, &grid_rows_, &grid_cols_, &coord_[0], &coord_[1]);
}

void BlockCyclicDistributionInformation::release_blacs_grid()
{
    if (owns_blacs_context_ && blacs_context_ >= 0)
    {
        Cblacs_gridexit(blacs_context_);
        blacs_context_ = kUnownedBlacsContext;
        owns_blacs_context_ = false;
    }
}

void BlockCyclicDistributionInformation::set_distribution_info(const int global_rows,
                                                               const int global_cols,
                                                               const int block_size)
{
    is_serial_ = false;
    block_size_ = block_size;

    // number of local rows and columns
    const int zero = 0;
    number_of_local_rows_ = numroc_(&global_rows, &block_size, &coord_[0], &zero, &grid_rows_);
    number_of_local_cols_ = numroc_(&global_cols, &block_size, &coord_[1], &zero, &grid_cols_);
    number_of_local_elements_ = static_cast<int64_t>(number_of_local_rows_) * number_of_local_cols_;

    // initialize the ScaLAPACK descriptor
    int info = 0;
    int lld = std::max(number_of_local_rows_, 1);
    const int descriptor_block_size = block_size_;
    const int descriptor_blacs_context = blacs_context_;
    descinit_(descriptor_.data(),
              &global_rows,
              &global_cols,
              &descriptor_block_size,
              &descriptor_block_size,
              &zero,
              &zero,
              &descriptor_blacs_context,
              &lld,
              &info);

    // generate the global-to-local and local-to-global index maps
    local_to_global_row_.resize(number_of_local_rows_);
    global_to_local_row_ = std::vector<int>(global_rows, -1);
    for (int i = 0; i < number_of_local_rows_; ++i)
    {
        local_to_global_row_[i]
            = (i / block_size * grid_rows_ + coord_[0]) * block_size + i % block_size;
        global_to_local_row_[local_to_global_row_[i]] = i;
    }

    local_to_global_col_.resize(number_of_local_cols_);
    global_to_local_col_ = std::vector<int>(global_cols, -1);
    for (int j = 0; j < number_of_local_cols_; ++j)
    {
        local_to_global_col_[j]
            = (j / block_size * grid_cols_ + coord_[1]) * block_size + j % block_size;
        global_to_local_col_[local_to_global_col_[j]] = j;
    }
}

int BlockCyclicDistributionInformation::initialize(const int global_rows,
                                                   const int global_cols,
                                                   const int block_size,
                                                   const MPI_Comm comm)
{
    return initialize(global_rows, global_cols, block_size, comm, false);
}

int BlockCyclicDistributionInformation::initialize(const int global_rows,
                                                   const int global_cols,
                                                   const int block_size,
                                                   const MPI_Comm comm,
                                                   const bool swap_grid)
{
    validate_global_shape(global_rows, global_cols);
    validate_block_size(block_size);

    release_blacs_grid();
    initialize_blacs_grid(comm, swap_grid);
    owns_blacs_context_ = true;
    set_distribution_info(global_rows, global_cols, block_size);
    return number_of_local_rows_ == 0 || number_of_local_cols_ == 0;
}

int BlockCyclicDistributionInformation::initialize(const int global_rows,
                                                   const int global_cols,
                                                   const int number_of_row_processors,
                                                   const int number_of_column_processors,
                                                   const int block_size,
                                                   const MPI_Comm communicator)
{
    validate_global_shape(global_rows, global_cols);
    validate_block_size(block_size);
    if (number_of_row_processors <= 0)
    {
        throw std::invalid_argument("number_of_row_processors must be positive");
    }
    if (number_of_column_processors <= 0)
    {
        throw std::invalid_argument("number_of_column_processors must be positive");
    }

    const int communicator_size = Parallel_Common::comm_size(communicator);
    const long long grid_size
        = static_cast<long long>(number_of_row_processors) * number_of_column_processors;
    if (grid_size != communicator_size)
    {
        throw std::invalid_argument(
            "number_of_row_processors * number_of_column_processors must equal communicator size");
    }

    release_blacs_grid();
    grid_rows_ = number_of_row_processors;
    grid_cols_ = number_of_column_processors;
    blacs_context_ = Csys2blacs_handle(communicator);
    char order = 'R';
    Cblacs_gridinit(&blacs_context_, &order, grid_rows_, grid_cols_);
    Cblacs_gridinfo(blacs_context_, &grid_rows_, &grid_cols_, &coord_[0], &coord_[1]);
    owns_blacs_context_ = true;
    set_distribution_info(global_rows, global_cols, block_size);
    return number_of_local_rows_ == 0 || number_of_local_cols_ == 0;
}

int BlockCyclicDistributionInformation::bind_descriptor(const int global_rows,
                                                        const int global_cols,
                                                        const int block_size,
                                                        const int blacs_context)
{
    validate_global_shape(global_rows, global_cols);
    validate_block_size(block_size);
    if (blacs_context < 0)
    {
        throw std::invalid_argument("blacs_context must be valid");
    }

    // Reusing the context owned by this object must not destroy the grid or
    // drop ownership; otherwise the grid about to be reused is released.
    if (blacs_context == blacs_context_ && owns_blacs_context_)
    {
        Cblacs_gridinfo(blacs_context, &grid_rows_, &grid_cols_, &coord_[0], &coord_[1]);
        set_distribution_info(global_rows, global_cols, block_size);
        return number_of_local_rows_ == 0 || number_of_local_cols_ == 0;
    }
    release_blacs_grid();
    blacs_context_ = blacs_context;
    owns_blacs_context_ = false;
    Cblacs_gridinfo(blacs_context, &grid_rows_, &grid_cols_, &coord_[0], &coord_[1]);
    set_distribution_info(global_rows, global_cols, block_size);
    return number_of_local_rows_ == 0 || number_of_local_cols_ == 0;
}

int BlockCyclicDistributionInformation::init(const int global_rows,
                                             const int global_cols,
                                             const int block_size,
                                             const MPI_Comm comm,
                                             const bool swap_grid)
{
    return initialize(global_rows, global_cols, block_size, comm, swap_grid);
}

int BlockCyclicDistributionInformation::set(const int global_rows,
                                            const int global_cols,
                                            const int block_size,
                                            const int blacs_context)
{
    return bind_descriptor(global_rows, global_cols, block_size, blacs_context);
}
#endif

void BlockCyclicDistributionInformation::initialize_serial(const int global_rows, const int global_cols)
{
    validate_global_shape(global_rows, global_cols);

    block_size_ = 1;
    grid_rows_ = 1;
    grid_cols_ = 1;
    coord_[0] = 0;
    coord_[1] = 0;
    number_of_local_rows_ = global_rows;
    number_of_local_cols_ = global_cols;
    number_of_local_elements_
        = static_cast<int64_t>(number_of_local_rows_) * number_of_local_cols_;
    local_to_global_row_.resize(number_of_local_rows_);
    local_to_global_col_.resize(number_of_local_cols_);
    std::iota(local_to_global_row_.begin(), local_to_global_row_.end(), 0);
    std::iota(local_to_global_col_.begin(), local_to_global_col_.end(), 0);
    global_to_local_row_ = local_to_global_row_;
    global_to_local_col_ = local_to_global_col_;
    is_serial_ = true;
#ifdef __MPI
    release_blacs_grid();
    // A serial layout must not reference any BLACS grid: for borrowers
    // release_blacs_grid() is a no-op, so clear the handle explicitly;
    // otherwise the object keeps a reference to the old grid, which also
    // dangles once the owner destroys it.
    blacs_context_ = kUnownedBlacsContext;
#endif
}

#ifdef __MPI
std::shared_ptr<const BlockCyclicDistributionInformation> make_frozen_distribution(
    const int global_rows,
    const int global_cols,
    const int block_size,
    const MPI_Comm communicator)
{
    validate_global_shape(global_rows, global_cols);
    validate_block_size(block_size);

    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->initialize(global_rows, global_cols, block_size, communicator);
    return distribution;
}

std::shared_ptr<const BlockCyclicDistributionInformation> make_frozen_distribution(
    const int global_rows,
    const int global_cols,
    const int number_of_row_processors,
    const int number_of_column_processors,
    const int block_size,
    const MPI_Comm communicator)
{
    validate_global_shape(global_rows, global_cols);
    validate_block_size(block_size);

    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->initialize(global_rows,
                             global_cols,
                             number_of_row_processors,
                             number_of_column_processors,
                             block_size,
                             communicator);
    return distribution;
}

std::shared_ptr<const BlockCyclicDistributionInformation> make_frozen_distribution(
    const int global_rows,
    const int global_cols,
    const int block_size,
    const int blacs_context)
{
    validate_global_shape(global_rows, global_cols);
    validate_block_size(block_size);
    if (blacs_context < 0)
    {
        throw std::invalid_argument("blacs_context must be valid");
    }

    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->bind_descriptor(global_rows, global_cols, block_size, blacs_context);
    return distribution;
}
#endif

std::shared_ptr<const BlockCyclicDistributionInformation> make_frozen_distribution(
    const int global_rows,
    const int global_cols)
{
    validate_global_shape(global_rows, global_cols);

    std::shared_ptr<BlockCyclicDistributionInformation> distribution
        = std::make_shared<BlockCyclicDistributionInformation>();
    distribution->initialize_serial(global_rows, global_cols);
    return distribution;
}
