import qctoolkit as qtk
import os, re, shutil, copy
from qctoolkit.QM.pseudo.pseudo import PP

def Al1st(qminp, **setting):
  assert 'ref_dir' in setting
  assert os.path.exists(setting['ref_dir'])

  inpClass = qtk.QM.general_io.GenericQMInput
  is_qminp = issubclass(qminp.__class__, inpClass)
  if not is_qminp:
    if type(qminp) is str:
      qminp = qtk.Molecule(qminp, **setting)
    elif type(qminp) is qtk.Molecule:
      qminp = qtk.QMInp(qminp, **setting)
  if 'program' not in setting:
    setting['program'] = qtk.setting.qmcode
  qminp.setting['scf_step'] = 1

  name = qminp.molecule.name
  if 'out_dir' in setting:
    name = setting['out_dir']
    del setting['out_dir']

  if setting['program'] == 'cpmd':
    setting['restart'] = True
    rst = os.path.join(setting['ref_dir'], 'RESTART')
    assert os.path.exists(rst)
    if 'dependent_files' in setting:
      setting['dependent_files'].append(rst)
    else:
      setting['dependent_files'] = [rst]

  if setting['program'] == 'bigdft':
    pass

  if setting['program'] == 'nwchem':
    pass

  qmout = qminp.run(name, **setting)
  return qmout

def mutatePP(pp1, pp2, fraction):
  pp1 = pp1*(1-fraction)
  pp2 = pp2*fraction
  pp = pp1 + pp2
  if fraction > 0.5:
    pp.param['Z'] = pp2.param['Z']

  return pp
