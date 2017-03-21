import multiprocessing as mp
from psutil import virtual_memory
import os
import copy

__version__ = '0.0.13'
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
mpi_flags = []

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
               'link_dep': False,
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
dcscp_list = ['H', 'He', 'C', 'N', 'O', 'Ne', 'P', 'S', 'Kr']
libxc_dict = {
  'pbe': ['XC_GGA_X_PBE', 'XC_GGA_C_PBE'],
  'blyp': ['XC_GGA_X_B88', 'XC_GGA_C_LYP'],
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
#cpmd_pp_url = 'http://cp2k.web.psi.ch/potentials/Goedecker-Teter-Hutter/cpmd/'
cpmd_pp_url = 'https://sourceforge.net/p/cp2k/code/HEAD/tree/trunk/potentials/Goedecker/cpmd'
#cpmd_dcacp_url = 'http://lcbc.epfl.ch/files/content/sites/lcbc/files/DCACPs/download/SG/'
cpmd_dcacp_url = 'http://lcbc.epfl.ch/page-71135-en.html'
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
#bigdft_pp_url = 'http://cp2k.web.psi.ch/potentials/Goedecker-Teter-Hutter/abinit/'
bigdft_pp_nlcc_url = 'http://bigdft.org/Wiki/index.php?title=New_Soft-Accurate_NLCC_pseudopotentials'
bigdft_pp_url = 'https://sourceforge.net/p/cp2k/code/HEAD/tree/trunk/potentials/Goedecker/abinit/'
abinit_exe = 'abinit'
abinit_cut3d_exe = 'cut3d'
abinit_f2b_exe = 'fold2Bloch'
#bigdft_pp = '/home/samio/Works/PhD/packages/BigDFT/PP'
gaussian_exe = 'g09'
gaussian_cubegen_exe = 'cubegen'
gaussian_formchk_exe = 'formchk'


program_dict = {
  'cpmd': cpmd_exe,
  'vasp': vasp_exe,
  'espresso': espresso_exe + ' <',
  'abinit': abinit_exe + ' <',
  'bigdft': bigdft_exe,
  'nwchem': nwchem_exe,
  'gaussian': gaussian_exe,
}

popleList = [
  "3-21g",
  '3-21g*',
  '3-21++g',
  '3-21++g*',
  '6-31g',
  '6-31g*',
  '6-31g**',
  '6-31+g',
  '6-31+g*',
  '6-31+g**',
  '6-31++g',
  '6-31++g*',
  '6-31++g**',
  '6-311g',
  '6-311g*',
  '6-311g**',
  '6-311+g',
  '6-311+g*',
  '6-311+g**',
  '6-311++g',
  '6-311++g*',
  '6-311++g**',
]

ccList = [
  'cc-pVDZ',
  'cc-pVTZ',
  'cc-pVDZ',
  'cc-pVQZ',
  'cc-pV5Z',
  'cc-pV6Z',
  'cc-pV8Z',
  'cc-pV9Z',
  'aug-cc-pVDZ',
  'aug-cc-pVTZ',
  'aug-cc-pVDZ',
  'aug-cc-pVQZ',
  'aug-cc-pV5Z',
  'aug-cc-pV6Z',
]

defList = [
  'Def2-ECP',
  'Def2-SV',
  'Def2-SVP',
  'Def2-SVPD',
  'Def2-TZVP',
  'Def2-TZVPD',
  'Def2-TZVPP',
  'Def2-TZVPPD',
  'Def2-QZVP',
  'Def2-QZVPD',
  'Def2-QZVPP',
  'Def2-QZVPPD',
]

basis = [popleList, ccList, defList]
