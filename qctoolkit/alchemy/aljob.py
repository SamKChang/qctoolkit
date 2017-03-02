import qctoolkit as qtk
import os, re, shutil, copy, glob
from qctoolkit.QM.pseudo.pseudo import PP
import universal as univ

def Al1st(qminp, runjob=False, **setting):
  assert 'ref_dir' in setting
  assert os.path.exists(setting['ref_dir'])
  setting['ref_dir'] = os.path.abspath(setting['ref_dir'])

  if 'runjob' not in setting:
    setting['runjob'] = False

  qminp = copy.deepcopy(univ.toInp(qminp, **setting))
  if hasattr(qminp.setting, 'save_restart'):
    qminp.setting['save_restart'] = False

  name = qminp.molecule.name
  if 'out_dir' in setting:
    name = setting['out_dir']
    del setting['out_dir']
    qminp.molecule.name = name
    qminp.setting['root_dir'] = name

  if qminp.setting['program'] == 'cpmd':
    setting['restart'] = True
    setting['scf_step'] = 1
    rst = os.path.join(setting['ref_dir'], 'RESTART')
    assert os.path.exists(rst)
    rst = os.path.abspath(rst)
    if 'dependent_files' in setting:
      setting['dependent_files'].append(rst)
    else:
      setting['dependent_files'] = [rst]

  elif qminp.setting['program'] == 'abinit':
    setting['restart'] = True
    setting['scf_step'] = 0
    rstList = glob.glob(os.path.join(setting['ref_dir'], '*o_*WFK'))
    assert len(rstList) > 0
    rstSrc = rstList[-1]
    if 'dependent_file' in setting:
      setting['dependent_files'].append([rstSrc, name+'i_WFK'])
    else:
      setting['dependent_files'] = [[rstSrc, name+'i_WFK']]

    denList = glob.glob(os.path.join(setting['ref_dir'], '*o_*DEN'))
    if len(denList) > 0:
      setting['restart_density'] = True
      denSrc = denList[-1]
      setting['dependent_files'].append([denSrc, name+'i_DEN'])

  elif qminp.setting['program'] == 'espresso':
    wfn = glob.glob(setting['ref_dir'] + '/*.wfc[0-9]*')
    if 'threads' not in setting or setting['threads'] != len(wfn):
      qtk.warning('threads must be the same as ref_dir, reset to %d'\
                  % len(wfn))
      setting['threads'] = len(wfn)
    setting['restart'] = True
    setting['scf_step'] = 1
    rst = glob.glob(setting['ref_dir'] + '/*.restart*')
    save = glob.glob(setting['ref_dir'] + '/*.save')
    if 'dependent_files' in setting:
      for lst in [wfn, rst, save]:
        setting['dependent_files'].extend(lst)
    else:
      setting['dependent_files'] = wfn
      for lst in [rst, save]:
        setting['dependent_files'].extend(lst)

  elif qminp.setting['program'] == 'vasp':
    setting['restart'] = True
    setting['scf_step'] = 1
    wfn_list = glob.glob(setting['ref_dir'] + '/WAVECAR')
    assert len(wfn_list) == 1
    wfn = wfn_list[0]
    if 'dependent_files' in setting:
      setting['dependent_files'].append(wfn)
    else:
      setting['dependent_files'] = [wfn]

  elif qminp.setting['program'] == 'bigdft':
    pass

  elif qminp.setting['program'] == 'nwchem':
    pass

  if setting['runjob']:
    qmout = qminp.run(name, **setting)
    return qmout
  else:
    qminp.molecule.name = name
    qminp.setting.update(setting)
    new_inp = qtk.QMInp(qminp.molecule, **qminp.setting)
    return new_inp

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
