import qctoolkit as qtk
import os, re, shutil, copy, glob
from qctoolkit.QM.pseudo.pseudo import PP
import universal as univ

def Al1st(qminp, **setting):
  assert 'ref_dir' in setting
  assert os.path.exists(setting['ref_dir'])

  qminp = univ.toInp(qminp, **setting)

  name = qminp.molecule.name
  if 'out_dir' in setting:
    name = setting['out_dir']
    del setting['out_dir']
    qminp = copy.deepcopy(qminp)
    qminp.molecule.name = name

  if qminp.setting['program'] == 'cpmd':
    setting['restart'] = True
    setting['scf_step'] = 1
    rst = os.path.join(setting['ref_dir'], 'RESTART')
    assert os.path.exists(rst)
    if 'dependent_files' in setting:
      setting['dependent_files'].append(rst)
    else:
      setting['dependent_files'] = [rst]

  elif qminp.setting['program'] == 'espresso':
    wfn = glob.glob(setting['ref_dir'] + '/*.wfc[0-9]*')
    if 'threads' not in setting or setting['threads'] != len(wfn):
      qtk.warning('threads must be the same as ref_dir, reset to %d'\
                  % len(wfn))
      setting['threads'] = len(wfn)
    setting['restart'] = True
    setting['scf_step'] = 1

  elif qminp.setting['program'] == 'bigdft':
    pass

  elif qminp.setting['program'] == 'nwchem':
    pass

  #qminp.write(name, **setting)
  qmout = qminp.run(name, **setting)
  return qmout

def mutatePP(pp1, pp2, fraction):
  if type(pp1) is str:
    if pp1.upper() == 'VOID':
      pp1 = PP()
    else:
      pp1 = PP(pp1)
  if type(pp2) is str:
    if pp2.upper() == 'VOID':
      pp2 = PP()
    else:
      pp2 = PP(pp2)
  pp1 = pp1*(1-fraction)
  pp2 = pp2*fraction
  pp = pp1 + pp2
  if pp1.param['Z']*pp2.param['Z'] > 0:
    if fraction > 0.5:
      pp.param['Z'] = pp2.param['Z']
  else:
    if pp1.param['Z'] == 0:
      pp.param['Z'] = pp2.param['Z']
    else:
      pp.param['Z'] = pp1.param['Z']

  return pp
