import multiprocessing as mp
from psutil import virtual_memory

__version__ = '0.1'
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

# VASP setting
vasp_pp = '/home/samio/Works/PhD/packages/VASP/PP'
vasp_exe = 'vasp'

# QM executables
cpmd_exe = 'cpmd.x'
cpmd_cpmd2cube = 'cpmd2cube.x'
cpmd_pp = ''
