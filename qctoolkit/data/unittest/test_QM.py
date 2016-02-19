from setup_test import *
from nose.plugins.skip import Skip, SkipTest
import numpy as np

pw_list = ['cpmd', 'vasp']
g_list = ['nwchem']
modes = ['single_point', 'geopt']
pw_theory = ['pbe', 'blyp', 'lda']
g_theory = ['pbe', 'pbe0', 'blyp', 'b3lyp', 'bp91', 'bp86', 'pw91',
            'rhf', 'rohf', 'uhf', 
            'mp2', 'ccsd', 'ccsdt',
           ]
tce = [
        'mp2', 'mp3', 'mp4',
        'ccsd', 'ccsdt', 'lccsd',
        'cisd', 'cisdt',
      ]

tmp_str = 'qct_test_qm_'

def test_general_inp_render():
  mol = setup(mol='h2o.xyz')[0]
  tmp_name = tmp_str + 'inp_render_' + str(os.getpid())
  def testRun(codes, theories):
    for code in codes:
      for theory in theories:
        for mode in modes:
          inp = qtk.QMInp(mol, mode=mode, theory=theory, program=code)
          inp.write()
          tmp_obj, tmp_inp = inp.write(tmp_name)
          try:
            os.remove(tmp_inp)
          except OSError:
            shutil.rmtree(tmp_inp)
  testRun(pw_list, pw_theory)
  testRun(g_list, g_theory)

def test_h2_pbe_allcode():
  if qtk.setting.run_qmtest:
    codes = ['nwchem', 'cpmd', 'vasp']
    mol = setup(mol='h2.xyz')[0]
    for i in range(len(codes)):
      code = codes[i]
      name = tmp_str + code
      inp = qtk.QMInp(mol, program=code)
      out = inp.run(name, threads = 2)
  else:
    gn = '\033[92m'
    ec = '\033[0m'
    raise SkipTest("\n %sSkipping qm tests%s for speed! To turn on, " \
                   % (gn, ec) + "set setting.run_qmtest=True")

def test_h2_tce_singlet():
  if qtk.setting.run_qmtest:
    mol = setup(mol='h2.xyz')[0]
    Et = []
    for cc in tce:
      name = tmp_str + cc + "_singlet"
      inp = qtk.QMInp(mol, theory=cc, program='nwchem')
      out = inp.run(name)
      Et.append(out.Et)
      print Et
  else:
    raise SkipTest

def test_h2_tce_doublet():
  if qtk.setting.run_qmtest:
    mol = setup(mol='h2.xyz')[0]
    mol.setChargeMultiplicity(-1, 2)
    Et = []
    for cc in tce:
      name = tmp_str + cc + "doublet"
      inp = qtk.QMInp(mol, theory=cc, program='nwchem')
      out = inp.run(name)
      Et.append(out.Et)
      print Et
  else:
    raise SkipTest

def test_cleanup():
  tmp_files = glob.glob(tmp_str + '*')
  for tmp in tmp_files:
    try:
      os.remove(tmp)
    except OSError:
      shutil.rmtree(tmp)
