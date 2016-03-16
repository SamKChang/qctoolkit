import qctoolkit as qtk
import qmcode.cpmd as cpmd
import qmcode.espresso as espresso
import qmcode.vasp as vasp
import qmcode.nwchem as nwchem
import qmcode.gaussian as gaussian
import qmcode.bigdft as bigdft
import qmcode.gaussian as gaussian

def QMInp(molecule, **kwargs):

  if type(molecule) is str:
    molecule = qtk.Molecule(molecule, **kwargs)

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode

  if kwargs['program'].lower() == 'cpmd':
    kwargs['extension'] = 'inp'
    return cpmd.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'espresso':
    kwargs['extension'] = 'inp'
    return espresso.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'vasp':
    return vasp.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'nwchem':
    kwargs['extension'] = 'inp'
    return nwchem.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'gaussian':
    kwargs['extension'] = 'com'
    return gaussian.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'bigdft':
    kwargs['extension'] = 'yaml'
    return bigdft.inp(molecule, **kwargs)

def QMOut(out=None, **kwargs):

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode

  if kwargs['program'].lower() == 'cpmd':
    return cpmd.out(out)
  elif kwargs['program'].lower() == 'vasp':
    return vasp.out(out)
  elif kwargs['program'].lower() == 'espresso':
    return espresso.out(out)
  elif kwargs['program'].lower() == 'nwchem':
    return nwchem.out(out)
  elif kwargs['program'].lower() == 'gaussian':
    return gaussian.out(out)
  elif kwargs['program'].lower() == 'bigdft':
    return bigdft.out(out)
  else:
    qtk.exit("program: %s not reconized" % kwargs['program'])
