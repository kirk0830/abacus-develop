CHANGELOG = [
    'Jul 14, 2024: initial version',
    'Aug  8, 2024: Fix: correct the behavior of profile1d function in cube_manipulator.py #4892',
    'Nov 25, 2025: Refactor & Feature: add the repeat function to cube_manipulator.py #????'
]

DESCRIPTION = [
    'A tool to manipulate the Gaussian cube file.',
    'A useful functionality has been supported: duplicate the charge density of the primitive cell'
    ', then repeat it to the supercell, this would be useful for the calculation of the supercell'
]

PERIODIC_TABLE = [
    'X', # dummy element, to avoid the index starts from 0
    'H', 'He', 
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 
    'K', 'Ca', 
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
    'Rb', 'Sr', 
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 
    'Cs', 'Ba', 
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 
    'Fr', 'Ra', 
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]
import sys
import unittest
import argparse
from enum import Enum
from copy import deepcopy
from itertools import product as itprod
from typing import Dict, Optional, Tuple
from pathlib import Path

import numpy as np

class GaussianCubeFile:
    '''
    a class to read and write Gaussian cube file.
    '''
    def __init__(self, fn: str | Path):
        '''initialize the cube file.
        
        Parameters
        ----------
        fn : str | Path
            the path to the cube file.
        '''
        self.fn_ = Path(fn)

        # check its existence
        if not self.fn_.exists():
            raise FileNotFoundError(f'{self.fn_} does not exist.')
        
        # parse the file
        try:
            self.read(self.fn_)
        except Exception as e:
            raise RuntimeError(f'failed to read {self.fn_}') from e

        # report the file information
        info  = f'GaussianCubeFile: the file {self.fn_} has been parsed.\n'
        info += f'{self}\n'
        print(info, flush=True)

    def __str__(self) -> str:
        '''the user-readable information'''

        return '\n'.join([
             'PIXEL',
            f'Pixel dimension: {self.nx_}x{self.ny_}x{self.nz_}',
            f'Pixel: x=({",".join([f"{x:>12.4e}" for x in self.pxl_[0,:]])})',
            f'       y=({",".join([f"{x:>12.4e}" for x in self.pxl_[1,:]])})',
            f'       z=({",".join([f"{x:>12.4e}" for x in self.pxl_[2,:]])})',
             'in unit of Bohr',
             '',
             'CELL',
            f'Cell:  x=({",".join([f"{x:>12.4e}" for x in self.cell_[0,:]])})',
            f'       y=({",".join([f"{x:>12.4e}" for x in self.cell_[1,:]])})',
            f'       z=({",".join([f"{x:>12.4e}" for x in self.cell_[2,:]])})',
             'in unit of Bohr',
             '',
             'ATOM',
            f'Number of atoms: {self.nat_}',
             'Elem     Charge        X           Y           Z',
            '\n'.join([f'{PERIODIC_TABLE[an]:4} {ac:11.6f} {c[0]:11.6f} {c[1]:11.6f} {c[2]:11.6f}'
                       for an, ac, c in zip(self.atomic_number_, self.atomic_charge_, self.tau_)]),
             'in unit of Bohr',
        ])

    def write(self, fn: str | Path, **kwargs) -> None:
        '''write the cube file.
        
        Parameters
        ----------
        fn : str | Path
            the path to the output cube file
        kwargs : dict
            optional arguments for writing the cube file:
            - ndigits: int (default=6)
                number of decimal digits to write for float values

        '''
        nplace = kwargs.get('ndigits', 6)
        width = nplace + 7 # 6 digits + sign + decimal point

        def build_pixel(nx, ny, nz, pxl):
            '''pixel definition.'''
            assert isinstance(pxl, np.ndarray)
            assert pxl.shape == (3, 3)
            return f'{nx:4d} {pxl[0,0]:11.6f} {pxl[0,1]:11.6f} {pxl[0,2]:11.6f}\n' \
                   f'{ny:4d} {pxl[1,0]:11.6f} {pxl[1,1]:11.6f} {pxl[1,2]:11.6f}\n' \
                   f'{nz:4d} {pxl[2,0]:11.6f} {pxl[2,1]:11.6f} {pxl[2,2]:11.6f}\n'

        def build_atomic_data(atomic_number, atomic_charge, coord):
            '''atomic data'''
            assert len(atomic_number) == len(atomic_charge) == len(coord)
            assert isinstance(coord, np.ndarray)
            assert coord.shape == (len(atomic_number), 3)
            return '\n'.join([f'{an:4d} {ac:11.6f} {c[0]:11.6f} {c[1]:11.6f} {c[2]:11.6f}'
                              for an, ac, c in zip(atomic_number, atomic_charge, coord)])

        flattenrho = self.rho_.flatten()
        with open(fn, 'w') as f:
            f.write(f'{self.comment_[0]}\n')
            f.write(f'{self.comment_[1]}\n')
            f.write(f'{self.nat_:4d} {self.xo_:11.6f} {self.yo_:11.6f} {self.zo_:11.6f}\n')
            f.write(build_pixel(self.nx_, self.ny_, self.nz_, self.pxl_))
            f.write(build_atomic_data(self.atomic_number_, self.atomic_charge_, self.tau_))
            f.write('\n')
            for i in range(0, len(flattenrho), 6):
                f.write(' '.join([f'%{width}.{nplace}e' % x for x in flattenrho[i:i+6]]))
                f.write('\n')
        print('>>> Write a cube file:', fn, flush=True)

    def read(self, fn: str | Path):
        '''read the cube file.
        
        Parameters
        ----------
        fn : str | Path
            the path to the cube file to be read
            
        Returns
        -------
        None
        
        Raises
        ------
        FileNotFoundError
            if the specified file does not exist
        RuntimeError
            if there are any parsing errors while reading the file
            
        Notes
        -----
        The cube file format is parsed as follows:
        1. First 2 lines: comments
        2. Line 3: number of atoms and origin coordinates
        3. Lines 4-6: grid dimensions and voxel basis vectors
        4. Next N lines: atomic data (atomic number, charge, coordinates)
        5. Remaining lines: volumetric data (charge density values)
        '''
        # read all the contents into memory
        print('<<< Read a cube file:', fn, flush=True)
        ################################################################################
        #                                    READ
        ################################################################################
        with open(fn, 'r') as f:
            raw = f.readlines()

        ################################################################################
        #                               PARSE AND WASH
        ################################################################################
        # the first two lines are the plain text comments
        self.comment_ = [l.strip() for l in raw[:2]]

        # the total number of atoms, the origin that defines the coordinate system
        self.nat_ = int(raw[2].split()[0])
        self.xo_, self.yo_, self.zo_ = map(float, raw[2].split()[1:])

        # the pixel dimension and basis
        self.nx_ = int(raw[3].split()[0])
        self.ny_ = int(raw[4].split()[0])
        self.nz_ = int(raw[5].split()[0])
        # basis
        self.pxl_ = np.array([raw[i].split()[1:4] 
                              for i in range(3, 6)]).astype(float).reshape(3, 3)
        # the cell is defined based on the assumption that the pixel fills the cell
        self.cell_ = self.pxl_ * np.array([self.nx_, self.ny_, self.nz_])

        def parser(l: str) -> Tuple[int, float, list[float]]:
            '''parse the atomic data line.'''
            l = l.split()
            return int(l[0]), float(l[1]), list(map(float, l[2:]))
        self.atomic_number_, self.atomic_charge_, self.tau_ = \
            zip(*map(parser, raw[6:6+int(self.nat_)]))

        # the charge density
        self.rho_ = np.array([float(x) 
                              for l in raw[6+int(self.nat_):] 
                              for x in l.split()]).reshape(self.nx_, self.ny_, self.nz_)
        
        ################################################################################
        #                           CHARGE DENSITY CHECK
        ################################################################################
        print('\nAfter-read quick check: ')
        sumrho, vpxl = np.sum(self.rho_), abs(np.linalg.det(self.pxl_))
        if vpxl < 1e-10:
            raise RuntimeError(f'V(pixel) is too small: {vpxl:10.4e} Bohr^3')
        ne, ne0 = sumrho * vpxl, np.sum(self.atomic_charge_)
        print(f'      V(pixel): {vpxl:10.4e} Bohr^3')
        print(f'int(dr*rho(r)): {ne:10.4e} e\n')
        
        err_ne = min(abs(ne - ne0), abs(ne*2 - ne0))
        if err_ne > 1e-2:
            raise RuntimeError(f'Number of electrons inconsistent, error={err_ne}')

    def repeat(self, mx, my, mz) -> 'GaussianCubeFile':
        '''Repeat the cube file to create a supercell.

        Parameters
        ----------
        mx : int
            repeat times along x direction
        my : int
            repeat times along y direction
        mz : int
            repeat times along z direction

        Returns
        -------
        GaussianCubeFile
            a new cube file object representing the supercell

        Notes
        -----
        This method will:
        1. Multiply the grid dimensions by the repeat factors
        2. Replicate atomic charges and numbers accordingly
        3. Translate atomic coordinates to fill the supercell
        4. Tile the volumetric data to match the new dimensions
        '''

        print('Repeat the cube file to the supercell:', 
              f'({mx}, {my}, {mz})', flush=True)
        if any([x < 1 for x in (mx, my, mz)]):
            raise ValueError('Repeat factors must be positive integers')

        out = deepcopy(self)
        out.nx_ *= mx
        out.ny_ *= my
        out.nz_ *= mz
        out.atomic_charge_ = np.tile(out.atomic_charge_, mx * my * mz).tolist()
        out.atomic_number_ = np.tile(out.atomic_number_, mx * my * mz).tolist()
        out.nat_ *= mx * my * mz

        def repeat_tau(tau, mx, my, mz, cell):
            repeated = None
            for i, j, k in itprod(range(mx), range(my), range(mz)):
                tvec = np.sum(np.array([[i, j, k]]).T * cell, axis=0) # the translation vector
                print('Perform supercell translation on atoms:', 
                      f'({",".join([f"{x:>12.4e}" for x in tvec])})', 'R=', (i, j, k))
                taunew = tau + tvec
                repeated = np.vstack([repeated, taunew]) if repeated is not None else taunew
            # after the loop, we return
            return repeated.reshape(-1, 3)
        
        out.tau_ = repeat_tau(out.tau_, mx, my, mz, out.cell_)
        out.rho_ = np.tile(out.rho_, (mx, my, mz))
        return out

    def todict(self) -> Dict:
        '''convert the cube file to a dictionary.'''
        return {
            'comment': self.comment_,
            'natom': self.nat_,
            'origin': (self.xo_, self.yo_, self.zo_),
            'nx': self.nx_,
            'ny': self.ny_,
            'nz': self.nz_,
            'pxl': self.pxl_,
            'cell': self.cell_,
            'atomz': self.atomic_number_,
            'chg': self.atomic_charge_,
            'coords': self.tau_,
            'data': self.rho_
        }

    def __repr__(self) -> str:
        return f'<GaussianCubeFile: {self.fn_}>'

    def __eq__(self, other) -> bool:
        '''check if two cube files are equal.'''
        if not isinstance(other, GaussianCubeFile):
            return False
        return self.todict() == other.todict()

    def __ne__(self, other) -> bool:
        '''check if two cube files are not equal.'''
        if not isinstance(other, GaussianCubeFile):
            return True
        return not self == other

def calculate_1d_profile(cube: GaussianCubeFile, 
                         axis: int | str) -> np.ndarray:
    '''Calculate the 1D profile by integrating charge density along a specified axis.

    Parameters
    ----------
    cube : GaussianCubeFile
        The cube file object containing charge density data
    axis : int | str
        The axis along which to integrate ('x', 'y', 'z' or 0, 1, 2)

    Returns
    -------
    np.ndarray
        The integrated 1D profile array

    Notes
    -----
    The integration is performed by summing the charge density values
    along the specified axis, effectively collapsing the 3D grid into
    a 1D profile while preserving the total charge.
    '''

    axis = 'xyz'.index(axis) if isinstance(axis, str) else axis
    assert axis in range(3)
    rho = cube.rho_
    return np.sum(rho, axis=axis)

def calculate_2d_slice(cube: GaussianCubeFile, 
                       axis: int | str, 
                       taud: float) -> np.ndarray:
    '''perform the slice along a given axis.'''
    axis = 'xyz'.index(axis) if isinstance(axis, str) else axis
    assert axis in range(3)
    n = [cube.nx_, cube.ny_, cube.nz_]
    return cube.rho_.take(int(taud * n[axis]), axis=axis)

def calculate_axpy(c1: GaussianCubeFile, 
                   alpha: float,
                   c2: Optional[GaussianCubeFile] = None, 
                   beta: Optional[float] = None) -> np.ndarray:
    '''perform the axpy operation. Note: this function may raise the
    IndexError when the cube files have different pixel dimensions.'''
    rho = c1.rho_ * alpha
    if c2 is not None:
        rho += c2.rho_ * beta
    return rho

def entry() -> Dict[str, str]:
    '''
    parse the command line arguments. Supported options:
    
    Debug arguments:
    --unittest: run the unittest

    Compulsory arguments:
    -i / --input: the input cube file.
    -o / --output: the output file prefix

    Functional arguments:
    -r / --repeat: repeat the cube file, whose value will be like '2,3,4'
    -s / --slice: slice the cube file, whose value will be like 'x=0.5'
    -p / --profile: calculate the 1D profile, whose value will be like 'x'
    
    Special arguments for the functionality `axpy`:
    --plus: plus the cube file with the other, whose value should be the cube file path.
    --minus: minus the cube file with the other, whose value should be the cube file path.
    '''
    parser = argparse.ArgumentParser(
        description='\n'.join(DESCRIPTION),
        epilog='Changelog: ' + '|'.join(CHANGELOG)
    )

    #######################################################################################
    #                                   ARGUMENTS DEFINITION
    #######################################################################################
    parser.add_argument('--unittest', 
                        action='store_true', 
                        help='run the unittest')
    parser.add_argument('-i', 
                        '--input', 
                        type=str, 
                        default=None,
                        help='the input cube file.')
    parser.add_argument('-o', 
                        '--output', 
                        type=str, 
                        default=None, 
                        help='the output file prefix.')
    parser.add_argument('-r', 
                        '--repeat', 
                        type=str, 
                        default=None, 
                        help='repeat the cube file, whose value will be like \'2,3,4\'')
    parser.add_argument('-s', 
                        '--slice', 
                        type=str, 
                        default=None, 
                        help='slice the cube file, whose value will be like \'x=0.5\'')
    parser.add_argument('-p', 
                        '--profile', 
                        type=str, 
                        default=None, 
                        help='calculate the 1D profile, whose value will be like \'x\'')
    parser.add_argument('--plus', 
                        type=str, 
                        default=None, 
                        help='plus the cube file with the other, whose value should be the cube file path.')
    parser.add_argument('--minus', 
                        type=str, 
                        default=None, 
                        help='minus the cube file with the other, whose value should be the cube file path.')
    
    # return the dictionary rather than the namespace object
    return vars(parser.parse_args())

class WorkflowStatus(Enum):
    UNITTEST, NORMAL = range(2)

def main(plan: Dict[str, str | int | float]) -> WorkflowStatus:
    '''
    initialize the workflow, navigate to the proper sub-function.

    Parameters
    ----------
    plan : Dict[str, str | int | float]
        the plan dictionary

    Returns
    -------
    WorkflowStatus
        the workflow status
    '''
    if plan['unittest']:
        return WorkflowStatus.UNITTEST
    
    # we have to place the input and output check here
    if plan['input'] is None:
        raise ValueError('The input cube file is not specified.')
    if plan['output'] is None:
        raise ValueError('The output file prefix is not specified.')

    # read
    mycube = GaussianCubeFile(plan['input'])

    if plan['repeat'] is not None:
        mx, my, mz = map(int, plan['repeat'].split(','))
        repeated = mycube.repeat(mx, my, mz)
        repeated.write(f'{plan["output"]}-repeat.cube')

    if plan['slice'] is not None:
        axis, taud = plan['slice'].split('=')
        taud = float(taud)
        rho = calculate_2d_slice(mycube, axis, taud)
        np.savetxt(f'{plan["output"]}-2dslice.txt', rho)

    if plan['profile'] is not None:
        axis = plan['profile']
        rho = calculate_1d_profile(mycube, axis)
        np.savetxt(f'{plan["output"]}-1dprofile.txt', rho)

    if plan['plus'] is not None:
        othercube = GaussianCubeFile(plan['plus'])
        rho = calculate_axpy(mycube, 1.0, othercube, 1.0)
        newcube = deepcopy(mycube)
        newcube.rho_ = rho
        newcube.write(f'{plan["output"]}-sum.cube')

    if plan['minus'] is not None:
        othercube = GaussianCubeFile(plan['minus'])
        rho = calculate_axpy(mycube, 1.0, othercube, -1.0)
        newcube = deepcopy(mycube)
        newcube.rho_ = rho
        newcube.write(f'{plan["output"]}-diff.cube')

    return WorkflowStatus.NORMAL

class TestCubeManipulator(unittest.TestCase):

    here = Path(__file__).parent
    fcube = here.parent.parent / 'tests' / '01_PW' / '050_PW_CHG_mismatch' / 'chgs1.cube'

    def test_read(self):
        cube = GaussianCubeFile(self.fcube)
        self.assertIsInstance(cube, GaussianCubeFile)

if __name__ == '__main__':
    
    if main(entry()) == WorkflowStatus.UNITTEST:
        # clean all the arguments because the unittest framework will
        # parse the arguments again.
        sys.argv = sys.argv[:1]
        unittest.main(exit=True)