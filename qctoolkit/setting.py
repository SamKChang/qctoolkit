import multiprocessing as mp
from psutil import virtual_memory
import os

__version__ = '0.0.10'
quiet = False
debug = False
no_warning = False
cpu_count = mp.cpu_count()
memory = float(virtual_memory().total)/10**9
qmcode = 'cpmd'
mdcode = 'cpmd'
run_qmtest = False
download_pp = True
#run_qmtest = True

# geometry setting
bond_ratio = 1.1  # for Molecule.findBond function
box_margin = 3     # for planewave box setup

# MPI setting
mpistr = 'mpirun -np'
ompthreads = 1
ompstr = '-x OMP_NUM_THREADS=%d' % ompthreads

# QM executables
libgbasis = '/home/samio/src/science/nwchem-6.6/src/basis/libraries'
# default setup for qm jobs
qm_setup = {
             'threads': 1,
             'program' : qmcode,
             'theory' : 'pbe',
             'mode' : 'single_point',
             'scf_step' : 1000,
             'wf_convergence' : 1E-5,
             'geometry_convergence' : 1E-4,
             'fix_molecule' : True,
             'unit' : 'angstrom',
             'save_density' : False,
             'save_wf': False,
             'corner_cube': False,
             'save_restart' : False,
             'restart' : False,
             'prefix' : None,
             'suffix' : None,
             'extension' : None,
             'periodic': False,
             'restart': False,
             'setting_backup': False,
             'molecule_backup': False,
             'overwrite': False,
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
               'output' : False,
               'dependent_files' : [],
               'finalized': False,
             }
planewave_setup = {
                    'periodic': True,
                    'symmetry': 'orthorhombic',
                    'cutoff': 'mt',
                  }

dft_list = ['pbe', 'pbe0', 'blyp', 'b3lyp', 'bp91', 'bp86', 'pw91']
cc_list = ['mp2', 'mp3', 'mp4', 'ccsd', 'ccsdt', 
           'lccsd', 'cisd', 'cisdt']
dcacp_dict = {
  'C': 'C6H6',
  'H': 'H2',
  'N': 'N2',
  'O': 'CO2',
  'S': 'CS2',
}
_setting = os.path.realpath(__file__)
_root = os.path.split(_setting)[0]
# VASP setting
vasp_pp = '/home/samio/Works/PhD/packages/VASP/PP'
vasp_exe = 'vasp'
# CPMD setting
cpmd_exe = 'cpmd.x'
cpmd_cpmd2cube_exe = 'cpmd2cube.x'
cpmd_pp = os.path.join(_root, 'data/PP/cpmd')
cpmd_pp_url = 'http://cp2k.web.psi.ch/potentials/Goedecker-Teter-Hutter/cpmd/'
cpmd_dcacp_url = 'http://lcbc.epfl.ch/files/content/sites/lcbc/files/DCACPs/download/SG/'
# espresso setting
espresso_exe = 'pw.x'
espresso_cpmd2upf_exe = 'cpmd2upf.x'
espresso_pp = os.path.join(_root, 'data/PP/espresso')
espresso_pp_url = 'http://www.quantum-espresso.org/wp-content/uploads/upf_files/'
# NWChem setting
nwchem_exe = 'nwchem'
nwchem_mov2asc_exe = 'mov2asc'
# BigDFT setting
bigdft_exe = 'bigdft'
bigdft_tool_exe = 'bigdft-tool'
bigdft_pp = os.path.join(_root, 'data/PP/bigdft')
bigdft_pp_url = 'http://cp2k.web.psi.ch/potentials/Goedecker-Teter-Hutter/abinit/'
bigdft_pp_nlcc_url = 'http://bigdft.org/Wiki/index.php?title=New_Soft-Accurate_NLCC_pseudopotentials'
#bigdft_pp_url = 'https://sourceforge.net/p/cp2k/code/HEAD/tree/trunk/potentials/Goedecker/abinit/'
#bigdft_pp = '/home/samio/Works/PhD/packages/BigDFT/PP'
gaussian_exe = 'g09'
cubegen_exe = 'cubegen'

program_dict = {
  'cpmd': cpmd_exe,
  'vasp': vasp_exe,
  'espresso': espresso_exe + ' <',
  'bigdft': bigdft_exe,
  'nwchem': nwchem_exe,
  'gaussian': gaussian_exe,
}
