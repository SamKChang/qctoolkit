import qctoolkit as qtk
import qmcode.cpmd as cpmd
import qmcode.vasp as vasp

def QMInp(molecule, **kwargs):

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode

  if kwargs['program'].lower() == 'cpmd':
    return cpmd.inp(molecule, **kwargs)
  elif kwargs['program'].lower() == 'vasp':
    return vasp.inp(molecule, **kwargs)

def QMOut(out, **kwargs):

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode

  if kwargs['program'].lower() == 'cpmd':
    return cpmd.out(out)
  elif kwargs['program'].lower() == 'vasp':
    return vasp.out(out)
