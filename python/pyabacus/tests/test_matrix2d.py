from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

from pyabacus.ModuleBase import Matrix2dComplex128, Matrix2dFloat64
from pyabacus.ModuleBase._base_pack import has_mpi_support


if has_mpi_support:
    from mpi4py import MPI

    COMMUNICATOR = MPI.COMM_WORLD
    WORLD_RANK = COMMUNICATOR.Get_rank()
    WORLD_SIZE = COMMUNICATOR.Get_size()
else:
    COMMUNICATOR = None
    WORLD_RANK = 0
    WORLD_SIZE = 1


def _make_values(rows: int, cols: int, is_complex: bool) -> np.ndarray:
    real_values = np.arange(rows * cols, dtype=np.float64).reshape(rows, cols)
    if not is_complex:
        return real_values
    imaginary_values = (0.25 + real_values) * -0.5
    return real_values + 1j * imaginary_values


@pytest.mark.parametrize(
    ("matrix_type", "is_complex", "expected_dtype"),
    [
        (Matrix2dFloat64, False, np.dtype(np.float64)),
        (Matrix2dComplex128, True, np.dtype(np.complex128)),
    ],
)
def test_from_numpy_to_numpy_and_local_to_numpy(matrix_type, is_complex, expected_dtype):
    values = _make_values(4, 5, is_complex)
    matrix = matrix_type.from_numpy(values)

    assert isinstance(matrix, matrix_type)
    assert matrix.shape == (4, 5)
    assert matrix.global_shape == (4, 5)
    assert matrix.local_shape == (4, 5)
    assert not matrix.is_distributed

    exported = matrix.to_numpy()
    local = matrix.local_to_numpy()
    assert exported.dtype == expected_dtype
    assert local.dtype == expected_dtype
    assert exported.shape == (4, 5)
    assert local.shape == (4, 5)
    assert exported.flags.c_contiguous
    assert local.flags.c_contiguous
    np.testing.assert_array_equal(exported, values)
    np.testing.assert_array_equal(local, values)

    np.testing.assert_array_equal(
        matrix_type.from_numpy(np.asfortranarray(values)).to_numpy(),
        values,
    )


@pytest.mark.parametrize("matrix_type", [Matrix2dFloat64, Matrix2dComplex128])
def test_from_numpy_rejects_invalid_shapes(matrix_type):
    with pytest.raises(ValueError, match="2-dimensional"):
        matrix_type.from_numpy(np.arange(6, dtype=np.float64))

    with pytest.raises(ValueError, match="positive matrix dimensions"):
        matrix_type.from_numpy(np.empty((0, 3), dtype=np.float64))


def test_to_npy_writes_real_c_order_float64_file(tmp_path):
    if WORLD_SIZE != 1:
        pytest.skip("Host matrix file export is a single-process level test")

    values = _make_values(4, 5, False)
    output_path = tmp_path / "matrix_float64.npy"
    matrix = Matrix2dFloat64.from_numpy(values)

    assert matrix.to_npy(str(output_path), 0)
    loaded = np.load(output_path)
    assert loaded.dtype == np.dtype(np.float64)
    assert loaded.shape == (4, 5)
    assert loaded.flags.c_contiguous
    np.testing.assert_array_equal(loaded, values)


def test_to_npy_rejects_non_temporary_path():
    matrix = Matrix2dFloat64.from_numpy(np.ones((2, 2), dtype=np.float64))
    repository_output = Path(__file__).resolve().parents[3] / "forbidden_matrix.npy"

    with pytest.raises(ValueError, match="temporary directory"):
        matrix.to_npy(str(repository_output), 0)
    assert not repository_output.exists()


@pytest.mark.parametrize(
    ("matrix_type", "is_complex"),
    [(Matrix2dFloat64, False), (Matrix2dComplex128, True)],
)
def test_distributed_roundtrip_and_root_only_npy(matrix_type, is_complex, monkeypatch):
    if not has_mpi_support:
        pytest.skip("MPI bindings were not compiled in this build")

    shared_directory = None
    if WORLD_RANK == 0:
        shared_directory = tempfile.mkdtemp(prefix="pyabacus_matrix2d_")
    shared_directory = COMMUNICATOR.bcast(shared_directory, root=0)
    shared_path = Path(shared_directory)
    output_path = shared_path / "distributed_matrix.npy"
    root_rank = WORLD_SIZE - 1

    try:
        values = _make_values(9, 7, is_complex)
        grid_rows = 2 if WORLD_SIZE == 4 else 1
        grid_cols = WORLD_SIZE // grid_rows
        distributed = matrix_type(9, 7, grid_rows, grid_cols, 2)
        distributed.scatter(matrix_type.from_numpy(values))

        assert distributed.is_distributed
        assert distributed.global_shape == (9, 7)
        local = distributed.local_to_numpy()
        assert local.shape == distributed.local_shape
        assert local.flags.c_contiguous
        np.testing.assert_array_equal(distributed.to_numpy(), values)

        original_save = np.save

        def recording_save(file_path, array, *args, **kwargs):
            marker_path = shared_path / f"numpy_save_called_by_{WORLD_RANK}"
            marker_path.write_text("called", encoding="utf-8")
            return original_save(file_path, array, *args, **kwargs)

        monkeypatch.setattr(np, "save", recording_save)
        wrote_file = distributed.to_npy(str(output_path), root_rank)
        COMMUNICATOR.Barrier()

        assert wrote_file is (WORLD_RANK == root_rank)
        for source_rank in range(WORLD_SIZE):
            marker_path = shared_path / f"numpy_save_called_by_{source_rank}"
            assert marker_path.exists() is (source_rank == root_rank)

        loaded = np.load(output_path)
        assert loaded.shape == (9, 7)
        assert loaded.flags.c_contiguous
        np.testing.assert_array_equal(loaded, values)
    finally:
        COMMUNICATOR.Barrier()
        if WORLD_RANK == 0:
            shutil.rmtree(shared_directory)
        COMMUNICATOR.Barrier()
