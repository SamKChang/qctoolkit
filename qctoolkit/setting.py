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

# geometry setting
bond_ratio = 1.1
pw_margin = 3

# MPI setting
mpistr = 'mpirun -np'
ompthreads = 1
ompstr = '-x OMP_NUM_THREADS=%d' % ompthreads

# QM executables
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
               'restart_file' : None,
             }
# VASP setting
vasp_pp = '/home/samio/Works/PhD/packages/VASP/PP'
vasp_exe = 'vasp'
# CPMD setting
cpmd_exe = 'cpmd.x'
cpmd_cpmd2cube = 'cpmd2cube.x'
cpmd_pp = ''
# NWChem setting
nwchem_exe = 'nwchem'
nwchem_mov2asc = 'mov2asc'
# BigDFT setting
bigdft_exp = 'bigdft'
