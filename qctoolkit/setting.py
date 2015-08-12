import multiprocessing as mp
from psutil import virtual_memory

__version__ = '0.1'
quiet = False
no_warning = False
cpu_count = mp.cpu_count()
memory = float(virtual_memory().total)/10**9
