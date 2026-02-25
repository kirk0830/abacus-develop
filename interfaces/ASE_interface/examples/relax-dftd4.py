'''we investigate the influence of dispersion precision on predicting the 
interlayer distance of graphene

to make this example executable, please install the dftd4 package by:
conda install -c conda-forge dftd4
conda install -c conda-forge dftd4-python
'''
import tempfile
from pathlib import Path # a more Pythonic alternative to the os.path
here = Path(__file__).parent
# to the directory where the pseudopotential and orbital files are stored
# In your case you change to the appropriate one
pporb = here.parent.parent.parent / 'tests' / 'PP_ORB'

from dftd4.ase import DFTD4
from abacuslite import Abacus, AbacusProfile
from ase.calculators.mixing import SumCalculator
from ase.atoms import Atoms
from ase.geometry import cellpar_to_cell
from ase.optimize import BFGS
from ase.constraints import FixCartesian # fix the X and Y

import numpy as np

bilayer = Atoms(
    symbols = ['C' for _ in range(4)],
    scaled_positions=np.array([[  0,   0, 1/3],
                               [1/3, 2/3, 1/3],
                               [  0,   0, 2/3],
                               [1/3, 2/3, 2/3]]),
    cell=cellpar_to_cell([2.47, 2.47, 8, 90, 90, 120]),
    pbc=True
)

aprof = AbacusProfile(
    command='mpirun -np 8 abacus',
    pseudo_dir=pporb,
    orbital_dir=pporb,
    omp_num_threads=1,
)

# calculate the interlayer distance at PBE level
def calculate_interlayer_distance(atoms: Atoms,
                                  calculator: Abacus) -> float:
    atoms = atoms.copy()
    atoms.set_constraint([FixCartesian(a=list(range(4)), mask=(True, True, False))])
    atoms.calc = calculator
    dyn = BFGS(atoms, logfile='-')
    dyn.run(fmax=0.01)
    # after relax, get the distance
    pos = atoms.get_positions()
    zuniq = np.unique(pos[:, 2])
    return float(zuniq.max() - zuniq.min())

common = {
    'profile': aprof,
    'pseudopotentials': {'C': 'C_ONCV_PBE-1.0.upf'},
    'basissets': {'C': 'C_gga_8au_100Ry_2s2p1d.orb'},
    'inp': {
        'calculation': 'scf',
        'nspin': 1,
        'basis_type': 'lcao',
        'ks_solver': 'genelpa',
        'ecutwfc': 100,
        'symmetry': 1,
        'dft_functional': 'pbe'
    },
    'kpts': {
        'mode': 'mp-sampling',
        'gamma-centered': True,
        'nk': (8, 8, 1),
        'kshift': (0, 0, 0)
    }
}

dist = {}
with tempfile.TemporaryDirectory() as jobdir:
    abacus = Abacus(
        directory=jobdir,
        **common
    )
    # equip with dftd4
    abacus = SumCalculator([abacus, DFTD4(method='PBE')])
    dist['dftd4'] = calculate_interlayer_distance(bilayer, abacus)

with tempfile.TemporaryDirectory() as jobdir:
    abacus = Abacus(
        directory=jobdir,
        **common
    )
    dist['none'] = calculate_interlayer_distance(bilayer, abacus)

with tempfile.TemporaryDirectory() as jobdir:
    dftd3 = common.copy()
    dftd3['inp']['vdw_method'] = 'd3_bj'
    abacus = Abacus(
        directory=jobdir,
        **dftd3
    )
    dist['dftd3'] = calculate_interlayer_distance(bilayer, abacus)

print(dist)