from setup_test import *
from nose.plugins.skip import Skip, SkipTest
import numpy as np
qtk.setting.no_warning = True
qtk.setting.quiet = True

pw_list = ['cpmd', 'vasp']
g_list = ['nwchem']
wl_list = ['bigdft']
modes = ['single_point', 'geopt']
pw_theory = ['pbe', 'blyp', 'lda']
g_theory = ['pbe', 'pbe0', 'blyp', 'b3lyp', 'bp91', 'bp86', 'pw91',
            'rhf', 'rohf', 'uhf', 
            'mp2', 'ccsd', 'ccsdt',
           ]
wl_theory = ['pbe', 'pbe0', 'pbesol']
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
          print '\n... doing %s %s %s' % (code, theory, mode)
          inp = qtk.QMInp(mol, mode=mode, theory=theory, program=code)
          inp.write()
          inp.write(tmp_name)
          shutil.rmtree(tmp_name)
  testRun(pw_list, pw_theory)
  testRun(g_list, g_theory)
  testRun(wl_list, wl_theory)

def test_general_inp_render_file():
  path = os.path.realpath(__file__)
  path = re.sub('[a-zA-Z0-9\._\-]*$', '', path)
  mol = glob.glob(path + 'test_data/molecules/h2o.xyz')[0]
  tmp_name = tmp_str + 'inp_frender_' + str(os.getpid())
  def testRun(codes, theories):
    for code in codes:
      for theory in theories:
        for mode in modes:
          print '\n... doing %s %s %s' % (code, theory, mode)
          inp = qtk.QMInp(mol, mode=mode, theory=theory, program=code)
          assert inp.setting['theory'] == theory
          assert inp.setting['program'] == code
          inp.write()
          inp.write(tmp_name)
          shutil.rmtree(tmp_name)
  testRun(pw_list, pw_theory)
  testRun(g_list, g_theory)
  testRun(wl_list, wl_theory)

def test_h2_pbe_allcode():
  if qtk.setting.run_qmtest:
    codes = ['nwchem', 'cpmd', 'vasp', 'bigdft']
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

def test_QMOut_operations():
  path = os.path.realpath(__file__)
  path = re.sub('[a-zA-Z0-9\._\-]*$', '', path)
  qmdict = {'nwchem':'n', 'cpmd':'c', 'vasp':'v', 'bigdft':'b'}
  unit_dict = {'Hartree':1, 'eV':27.211396132, 'kcal/mol':627.509469}
  root = os.path.join(path, 'test_data/qmout/')
  out = []
  Et = []
  for code, symb in qmdict.iteritems():
    stem = 'h2' + symb
    out_file = stem + '.out'
    path_dir = os.path.join(root, stem)
    path_out = os.path.join(path_dir, out_file)
    qmout = qtk.QMOut(path_out, program=code)
    out.append(qmout)
    Et.append(qmout.Et)

  new_list1 = []
  new_list2 = []
  new_E1 = []
  new_E2 = []
  for i in range(len(out)-1):
    new_list1.append(out[i+1] + out[i])
    new_list2.append(out[i+1] - out[i])
    new_E1.append(Et[i+1] + Et[i])
    new_E2.append(Et[i+1] - Et[i])

  for i in range(len(new_list1)):
    out1 = new_list1[i]
    out2 = new_list2[i]
    E1 = new_E1[i]
    E2 = new_E2[i]
    for unit, factor in unit_dict.iteritems():
      out1.inUnit(unit)
      out1.inUnit(unit)
      out2.inUnit(unit)
      out2.inUnit(unit)
      assert out1.Et == E1 * factor
      assert out2.Et == E2 * factor

def test_cleanup():
  tmp_files = glob.glob(tmp_str + '*')
  for tmp in tmp_files:
    try:
      os.remove(tmp)
    except OSError:
      shutil.rmtree(tmp)
