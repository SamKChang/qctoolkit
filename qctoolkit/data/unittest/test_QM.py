from setup_test import *

pw_list = ['cpmd', 'vasp']
g_list = ['nwchem']
modes = ['single_point', 'geopt']
pw_theory = ['pbe', 'blyp', 'lda']
g_theory = ['pbe', 'pbe0', 'blyp', 'b3lyp', 'bp91', 'bp86', 'pw91',
            'rhf', 'rohf', 'uhf', 
            'mp2', 'ccsd', 'ccsd(t)',
           ]

def test_general_io():
  mol = setup(mol='h2o.xyz')[0]
  def testRun(codes, theories):
    for code in codes:
      for theory in theories:
        for mode in modes:
          inp = qtk.QMInp(mol, mode=mode, theory=theory, program=code)
          inp.write()
  testRun(pw_list, pw_theory)
  testRun(g_list, g_theory)
