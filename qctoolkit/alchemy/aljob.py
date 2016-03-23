import qctoolkit as qtk
import os, re, shutil, copy, glob
from qctoolkit.QM.pseudo.pseudo import PP
import universal as univ

def Al1st(qminp, **setting):
  assert 'ref_dir' in setting
  assert os.path.exists(setting['ref_dir'])

  qminp = univ.toInp(qminp, **setting)
  qminp.setting['scf_step'] = 1

  name = qminp.molecule.name
  if 'out_dir' in setting:
    name = setting['out_dir']
    del setting['out_dir']

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
    setting['restart'] = True
    setting['scf_step'] = 1
    save_file = os.path.join(setting['ref_dir'], '*.save')
    save = glob.glob(save_file)
    print 'yo'
    print save
    try:
      save = save[0]
    except:
      qtk.exit('*.save folder not found.')
    save_folder = os.path.split(save)[1]
    print save_folder
    prefix = os.path.splitext(save_folder)[0]
    rst_files = os.path.join(setting['ref_dir'], prefix + '*')
    rst = glob.glob(rst_files)
    print rst
    if 'dependent_files' in setting:
      setting['dependent_files'].extend(rst)
    else:
      setting['dependent_files'] = rst

    print setting['dependent_files']
    # need to change pseudopotential name in pwscf.save

  elif qminp.setting['program'] == 'bigdft':
    pass

  elif qminp.setting['program'] == 'nwchem':
    pass

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
