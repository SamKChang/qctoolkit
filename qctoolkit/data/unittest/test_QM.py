from setup_test import *
from nose.plugins.skip import Skip, SkipTest
import numpy as np

pw_list = ['cpmd', 'vasp']
g_list = ['nwchem']
modes = ['single_point', 'geopt']
pw_theory = ['pbe', 'blyp', 'lda']
g_theory = ['pbe', 'pbe0', 'blyp', 'b3lyp', 'bp91', 'bp86', 'pw91',
            'rhf', 'rohf', 'uhf', 
            'mp2', 'ccsd', 'ccsd(t)',
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

def test_h2o_pbe():
  if qtk.setting.run_qmtest:
    codes = ['nwchem', 'cpmd', 'vasp']
    results = [-76.2982, -17.1334, -0.5220]
    mol = setup(mol='h2o.xyz')[0]
    for i in range(len(codes)):
      code = codes[i]
      Et = results[i]
      name = tmp_str + code
      inp = qtk.QMInp(mol, program=code)
      out = inp.run(name, threads = 2)
      print i
      print code
      print name
      print out.Et
      print Et
      print abs(out.Et - Et)
      assert abs(out.Et - Et) < 1E-4
  else:
    gn = '\033[92m'
    ec = '\033[0m'
    raise SkipTest("\n %sSkipping h2o_pbe%s for speed! To turn on, " \
                   % (gn, ec) + "set setting.run_qmtest=True")

def test_cleanup():
  tmp_files = glob.glob(tmp_str + '*')
  for tmp in tmp_files:
    try:
      os.remove(tmp)
    except OSError:
      shutil.rmtree(tmp)
