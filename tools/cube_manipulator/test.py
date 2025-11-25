#!/usr/bin/env python3
"""
Comprehensive unit tests for cube_manipulator.py

This test suite validates the functionality of the GaussianCubeFile class
and associated functions, including edge cases and error conditions.
"""

import unittest
import tempfile
import os
import numpy as np
from pathlib import Path
from copy import deepcopy
import sys

# Add the parent directory to the path to import the module
sys.path.insert(0, str(Path(__file__).parent))

from cube_manipulator import (
    GaussianCubeFile, 
    calculate_1d_profile, 
    calculate_2d_slice, 
    calculate_axpy,
    CHANGELOG,
    DESCRIPTION,
    PERIODIC_TABLE
)


class TestGaussianCubeFileBasic(unittest.TestCase):
    """Basic test cases for the GaussianCubeFile class."""

    def setUp(self):
        """Set up test fixtures before each test method."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_cube_file = self._create_simple_test_cube_file()

    def tearDown(self):
        """Clean up after each test method."""
        if os.path.exists(self.test_cube_file):
            os.remove(self.test_cube_file)
        os.rmdir(self.temp_dir)

    def _create_simple_test_cube_file(self) -> str:
        """Create a simple valid cube file for testing."""
        # Simple 2x2x2 cube with 1 atom, density sums to 1 electron
        content = """Simple test cube
Basic test file
1 0.0 0.0 0.0
2 1.0 0.0 0.0
2 0.0 1.0 0.0
2 0.0 0.0 1.0
1 1.0 0.0 0.0 0.0
0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125
"""
        file_path = os.path.join(self.temp_dir, 'simple_test.cube')
        with open(file_path, 'w') as f:
            f.write(content)
        return file_path

    def test_init_valid_file(self):
        """Test initialization with a valid cube file."""
        cube = GaussianCubeFile(self.test_cube_file)
        self.assertIsInstance(cube, GaussianCubeFile)
        self.assertEqual(cube.nat_, 1)
        self.assertEqual(cube.nx_, 2)
        self.assertEqual(cube.ny_, 2)
        self.assertEqual(cube.nz_, 2)
        self.assertEqual(len(cube.atomic_number_), 1)
        self.assertEqual(len(cube.atomic_charge_), 1)
        self.assertEqual(len(cube.tau_), 1)

    def test_init_nonexistent_file(self):
        """Test initialization with a non-existent file."""
        with self.assertRaises(FileNotFoundError):
            GaussianCubeFile('nonexistent.cube')

    def test_str_representation(self):
        """Test string representation of the cube file."""
        cube = GaussianCubeFile(self.test_cube_file)
        str_repr = str(cube)
        self.assertIn('PIXEL', str_repr)
        self.assertIn('CELL', str_repr)
        self.assertIn('ATOM', str_repr)
        self.assertIn('2x2x2', str_repr)

    def test_repr_representation(self):
        """Test repr representation of the cube file."""
        cube = GaussianCubeFile(self.test_cube_file)
        repr_str = repr(cube)
        self.assertIn('GaussianCubeFile', repr_str)

    def test_todict(self):
        """Test conversion to dictionary."""
        cube = GaussianCubeFile(self.test_cube_file)
        cube_dict = cube.todict()
        
        required_keys = ['comment', 'natom', 'origin', 'nx', 'ny', 'nz', 
                        'pxl', 'cell', 'atomz', 'chg', 'coords', 'data']
        for key in required_keys:
            self.assertIn(key, cube_dict)

    def test_repeat_valid(self):
        """Test valid repeat operation."""
        cube = GaussianCubeFile(self.test_cube_file)
        repeated = cube.repeat(2, 2, 2)
        
        # Check dimensions
        self.assertEqual(repeated.nx_, cube.nx_ * 2)
        self.assertEqual(repeated.ny_, cube.ny_ * 2)
        self.assertEqual(repeated.nz_, cube.nz_ * 2)
        self.assertEqual(repeated.nat_, cube.nat_ * 8)  # 2*2*2
        
        # Check atomic data
        self.assertEqual(len(repeated.atomic_number_), cube.nat_ * 8)
        self.assertEqual(len(repeated.atomic_charge_), cube.nat_ * 8)
        self.assertEqual(len(repeated.tau_), cube.nat_ * 8)

    def test_repeat_invalid_parameters(self):
        """Test repeat operation with invalid parameters."""
        cube = GaussianCubeFile(self.test_cube_file)
        
        # Test with zero
        with self.assertRaises(ValueError):
            cube.repeat(0, 2, 2)
        
        # Test with negative numbers
        with self.assertRaises(ValueError):
            cube.repeat(-1, 2, 2)

    def test_repeat_single(self):
        """Test repeat operation with (1,1,1)."""
        cube = GaussianCubeFile(self.test_cube_file)
        repeated = cube.repeat(1, 1, 1)
        
        # Should be identical to original
        self.assertEqual(repeated.nx_, cube.nx_)
        self.assertEqual(repeated.ny_, cube.ny_)
        self.assertEqual(repeated.nz_, cube.nz_)
        self.assertEqual(repeated.nat_, cube.nat_)


class TestCalculationFunctions(unittest.TestCase):
    """Test cases for calculation functions."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_cube_file = self._create_test_cube_file()
        self.cube = GaussianCubeFile(self.test_cube_file)

    def tearDown(self):
        """Clean up after tests."""
        if os.path.exists(self.test_cube_file):
            os.remove(self.test_cube_file)
        os.rmdir(self.temp_dir)

    def _create_test_cube_file(self) -> str:
        """Create a test cube file for calculation tests."""
        content = """Test cube for calculations
Generated for testing calculations
1 0.0 0.0 0.0
2 1.0 0.0 0.0
2 0.0 1.0 0.0
2 0.0 0.0 1.0
1 1.0 0.0 0.0 0.0
0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125
"""
        file_path = os.path.join(self.temp_dir, 'calc_test.cube')
        with open(file_path, 'w') as f:
            f.write(content)
        return file_path

    def test_calculate_1d_profile_string_axis(self):
        """Test 1D profile calculation with string axis."""
        profile_x = calculate_1d_profile(self.cube, 'x')
        profile_y = calculate_1d_profile(self.cube, 'y')
        profile_z = calculate_1d_profile(self.cube, 'z')
        
        # Check shapes - when summing along axis, we get 2D arrays
        self.assertEqual(profile_x.shape, (self.cube.ny_, self.cube.nz_))
        self.assertEqual(profile_y.shape, (self.cube.nx_, self.cube.nz_))
        self.assertEqual(profile_z.shape, (self.cube.nx_, self.cube.ny_))

    def test_calculate_1d_profile_int_axis(self):
        """Test 1D profile calculation with integer axis."""
        profile_0 = calculate_1d_profile(self.cube, 0)
        profile_1 = calculate_1d_profile(self.cube, 1)
        profile_2 = calculate_1d_profile(self.cube, 2)
        
        # Should be equivalent to string versions
        profile_x = calculate_1d_profile(self.cube, 'x')
        profile_y = calculate_1d_profile(self.cube, 'y')
        profile_z = calculate_1d_profile(self.cube, 'z')
        
        np.testing.assert_array_equal(profile_0, profile_x)
        np.testing.assert_array_equal(profile_1, profile_y)
        np.testing.assert_array_equal(profile_2, profile_z)

    def test_calculate_1d_profile_invalid_axis(self):
        """Test 1D profile calculation with invalid axis."""
        with self.assertRaises((AssertionError, ValueError)):
            calculate_1d_profile(self.cube, 'w')
        
        with self.assertRaises((AssertionError, ValueError)):
            calculate_1d_profile(self.cube, 3)

    def test_calculate_2d_slice_string_axis(self):
        """Test 2D slice calculation with string axis."""
        slice_x = calculate_2d_slice(self.cube, 'x', 0.5)
        slice_y = calculate_2d_slice(self.cube, 'y', 0.5)
        slice_z = calculate_2d_slice(self.cube, 'z', 0.5)
        
        # Check shapes
        self.assertEqual(slice_x.shape, (self.cube.ny_, self.cube.nz_))
        self.assertEqual(slice_y.shape, (self.cube.nx_, self.cube.nz_))
        self.assertEqual(slice_z.shape, (self.cube.nx_, self.cube.ny_))

    def test_calculate_2d_slice_int_axis(self):
        """Test 2D slice calculation with integer axis."""
        slice_0 = calculate_2d_slice(self.cube, 0, 0.5)
        slice_1 = calculate_2d_slice(self.cube, 1, 0.5)
        slice_2 = calculate_2d_slice(self.cube, 2, 0.5)
        
        # Should be equivalent to string versions
        slice_x = calculate_2d_slice(self.cube, 'x', 0.5)
        slice_y = calculate_2d_slice(self.cube, 'y', 0.5)
        slice_z = calculate_2d_slice(self.cube, 'z', 0.5)
        
        np.testing.assert_array_equal(slice_0, slice_x)
        np.testing.assert_array_equal(slice_1, slice_y)
        np.testing.assert_array_equal(slice_2, slice_z)

    def test_calculate_2d_slice_boundary_values(self):
        """Test 2D slice calculation with boundary values."""
        # Test with taud = 0 (first slice)
        slice_first = calculate_2d_slice(self.cube, 'x', 0.0)
        
        # Test with taud = 0.99 (last slice)
        slice_last = calculate_2d_slice(self.cube, 'x', 0.99)
        
        # Both should have correct shapes
        self.assertEqual(slice_first.shape, (self.cube.ny_, self.cube.nz_))
        self.assertEqual(slice_last.shape, (self.cube.ny_, self.cube.nz_))

    def test_calculate_2d_slice_invalid_axis(self):
        """Test 2D slice calculation with invalid axis."""
        with self.assertRaises((AssertionError, ValueError)):
            calculate_2d_slice(self.cube, 'w', 0.5)
        
        with self.assertRaises((AssertionError, ValueError)):
            calculate_2d_slice(self.cube, 3, 0.5)

    def test_calculate_axpy_single_cube(self):
        """Test axpy operation with a single cube."""
        result = calculate_axpy(self.cube, 2.0)
        expected = self.cube.rho_ * 2.0
        np.testing.assert_array_equal(result, expected)

    def test_calculate_axpy_two_cubes(self):
        """Test axpy operation with two cubes."""
        # Create a second cube with different data
        cube2 = deepcopy(self.cube)
        cube2.rho_ = cube2.rho_ + 1.0
        
        result = calculate_axpy(self.cube, 2.0, cube2, 3.0)
        expected = self.cube.rho_ * 2.0 + cube2.rho_ * 3.0
        np.testing.assert_array_equal(result, expected)

    def test_calculate_axpy_none_cube(self):
        """Test axpy operation with None cube."""
        result = calculate_axpy(self.cube, 2.0, None, None)
        expected = self.cube.rho_ * 2.0
        np.testing.assert_array_equal(result, expected)


class TestConstantsAndUtilities(unittest.TestCase):
    """Test cases for constants and utility functions."""

    def test_periodic_table_length(self):
        """Test that the periodic table has the expected length."""
        # Should have 119 elements (including dummy 'X')
        self.assertEqual(len(PERIODIC_TABLE), 119)

    def test_periodic_table_first_elements(self):
        """Test the first few elements of the periodic table."""
        expected_first = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne']
        self.assertEqual(PERIODIC_TABLE[:len(expected_first)], expected_first)

    def test_changelog_not_empty(self):
        """Test that CHANGELOG is not empty."""
        self.assertIsInstance(CHANGELOG, list)
        self.assertGreater(len(CHANGELOG), 0)

    def test_description_not_empty(self):
        """Test that DESCRIPTION is not empty."""
        self.assertIsInstance(DESCRIPTION, list)
        self.assertGreater(len(DESCRIPTION), 0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up after tests."""
        os.rmdir(self.temp_dir)

    def test_single_atom_cube(self):
        """Test cube file with single atom."""
        content = """Single atom cube
Test with one atom
1 0.0 0.0 0.0
2 1.0 0.0 0.0
2 0.0 1.0 0.0
2 0.0 0.0 1.0
1 1.0 0.0 0.0 0.0
0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125
"""
        file_path = os.path.join(self.temp_dir, 'single_atom.cube')
        with open(file_path, 'w') as f:
            f.write(content)
        
        try:
            cube = GaussianCubeFile(file_path)
            self.assertEqual(cube.nat_, 1)
            self.assertEqual(len(cube.atomic_number_), 1)
        finally:
            os.remove(file_path)


class TestErrorHandling(unittest.TestCase):
    """Test error handling and exception cases."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up after tests."""
        os.rmdir(self.temp_dir)

    def test_empty_cube_file(self):
        """Test handling of empty cube file."""
        empty_file = os.path.join(self.temp_dir, 'empty.cube')
        with open(empty_file, 'w') as f:
            f.write("")
        
        try:
            with self.assertRaises((IndexError, RuntimeError, ValueError)):
                GaussianCubeFile(empty_file)
        finally:
            os.remove(empty_file)

    def test_degenerate_cube_file(self):
        """Test cube file with degenerate pixel matrix."""
        content = """Degenerate cube file
Zero volume pixel
1 0.0 0.0 0.0
1 0.0 0.0 0.0
1 0.0 0.0 0.0
1 0.0 0.0 0.0
1 1.0 0.0 0.0 0.0
1.0
"""
        file_path = os.path.join(self.temp_dir, 'degenerate.cube')
        with open(file_path, 'w') as f:
            f.write(content)
        
        try:
            with self.assertRaises(RuntimeError):
                GaussianCubeFile(file_path)
        finally:
            os.remove(file_path)


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)