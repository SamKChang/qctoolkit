import multiprocessing as mp
from psutil import virtual_memory
import os

__version__ = '0.1.1'
quiet = False
no_warning = False
cpu_count = mp.cpu_count()
memory = float(virtual_memory().total)/10**9
qmcode = 'cpmd'
mdcode = 'cpmd'
run_qmtest = False
#run_qmtest = True

# geometry setting
bond_ratio = 1.1  # for Molecule.findBond function
pw_margin = 3     # for planewave box setup

# MPI setting
mpistr = 'mpirun -np'
ompthreads = 1
ompstr = '-x OMP_NUM_THREADS=%d' % ompthreads

# QM executables
libgbasis = '/home/samio/src/science/nwchem-6.6/src/basis/libraries'
# default setup for qm jobs
qm_setup = {
             'program' : qmcode,
             'theory' : 'pbe',
             'mode' : 'single_point',
             'scf_step' : 1000,
             'wf_convergence' : 1E-5,
             'geometry_convergence' : 1E-4,
             'fix_molecule' : True,
             'unit' : 'angstrom',
             'save_density' : False,
             'save_restart' : False,
             'restart' : False,
             'prefix' : None,
             'suffix' : None,
             'extension' : None,
             'periodic': False,
           }
# default setup for qmmm jobs
md_setup = {
             'T' : 298,
             'thermostat' : 'Nose-Hoover',
             'T_tolerance' : 0.1,
             'sample_period' : 10,
             'md_step' : 1000,
           }
# default setup for job output
file_setup = {
               'prefix' : None,
               'suffix' : None,
               'extension' : None,
               'root_dir' : None,
               'path' : os.getcwd(),
               'output' : False,
               'dependent_files' : [],
             }
planewave_setup = {
                    'periodic': True,
                    'symmetry': 'orthorhombic',
                    'cutoff': 'mt',
                  }

dft_list = ['pbe', 'pbe0', 'blyp', 'b3lyp', 'bp91', 'bp86', 'pw91']
cc_list = ['mp2', 'mp3', 'mp4', 'ccsd', 'ccsdt', 
           'lccsd', 'cisd', 'cisdt']
# VASP setting
vasp_pp = '/home/samio/Works/PhD/packages/VASP/PP'
vasp_exe = 'vasp'
# CPMD setting
cpmd_exe = 'cpmd.x'
cpmd_cpmd2cube = 'cpmd2cube.x'
cpmd_pp = '/home/samio/Works/PhD/packages/CPMD/PP'
# NWChem setting
nwchem_exe = 'nwchem'
nwchem_mov2asc = 'mov2asc'
# BigDFT setting
bigdft_exe = 'bigdft'
bigdft_pp = '/home/samio/Works/PhD/packages/BigDFT/PP'
