import qctoolkit as qtk
import qmcode.cpmd as cpmd
import qmcode.espresso as espresso
import qmcode.vasp as vasp
import qmcode.nwchem as nwchem
import qmcode.gaussian as gaussian
import qmcode.bigdft as bigdft
import qmcode.gaussian as gaussian
import sys, os

def QMInp(molecule, **kwargs):
  inp_dict = {
    'cpmd': [cpmd.inp, 'inp'],
    'vasp': [vasp.inp, ''],
    'espresso': [espresso.inp, 'inp'],
    'nwchem': [nwchem.inp, 'inp'],
    'gaussian': [gaussian.inp, 'com'],
    'bigdft': [bigdft.inp, 'ymal'],
  }

  if type(molecule) is str:
    molecule = qtk.Molecule(molecule, **kwargs)

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode
  p_str = kwargs['program'].lower()
  p = inp_dict[p_str][0]
  kwargs['extension'] = inp_dict[p_str][1]
  return p(molecule, **kwargs)

def QMOut(out=None, **kwargs):
  out_dict = {
    'cpmd': cpmd.out,
    'vasp': vasp.out,
    'espresso': espresso.out,
    'nwchem': nwchem.out,
    'gaussian': gaussian.out,
    'bigdft': bigdft.out,
  }

  if out is not None and not os.path.exists(out):
    qtk.warning("%s not found" % out)
    return qtk.QMOut()
  else:
    if 'program' in kwargs:
      p_str = kwargs['program']
      try:
        return out_dict[p_str](out)
      except:
        e = sys.exc_info()[0]
        qtk.warning("%s failed with message: %s" % (out, e))
        qout = qtk.QMOut()
        qout.path, qout.name = os.path.split(out)
        return qout
    else:
      for p in out_dict.itervalues():
        try:
          return p(out)
        except:
          pass
      qtk.warning("something wrong with output file, "+\
                  "pass 'program' eplicitly")
      
