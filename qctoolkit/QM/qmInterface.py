import qctoolkit as qtk
import qmcode.cpmd as cpmd
import qmcode.vasp as vasp
import qmcode.nwchem as nwchem
import qmcode.gaussian as gaussian

def QMInp(molecule, **kwargs):

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode

  if kwargs['program'].lower() == 'cpmd':
    return cpmd.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'vasp':
    return vasp.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'nwchem':
    return nwchem.inp(molecule, **kwargs)

def QMOut(out=None, **kwargs):

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode

  if kwargs['program'].lower() == 'cpmd':
    return cpmd.out(out)
  elif kwargs['program'].lower() == 'vasp':
    return vasp.out(out)
  elif kwargs['program'].lower() == 'nwchem':
    return nwchem.out(out)
  elif kwargs['program'].lower() == 'nwchem':
    return nwchem.out(out)
  elif kwargs['program'].lower() == 'gaussian':
    return gaussian.out(out)
  else:
    qtk.exit("program: %s not reconized" % kwargs['program'])
