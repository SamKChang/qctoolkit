import qctoolkit as qtk
import qmcode.cpmd as cpmd
import qmcode.espresso as espresso
import qmcode.vasp as vasp
import qmcode.abinit as abinit
import qmcode.nwchem as nwchem
import qmcode.gaussian as gaussian
import qmcode.bigdft as bigdft
import qmcode.gaussian as gaussian
import qmcode.hortonInterface as horton
import ofdft.main as ofdft
import sys, os, copy

def QMInp(molecule, **kwargs_in):
  kwargs = copy.deepcopy(kwargs_in)
  inp_dict = {
    'cpmd': [cpmd.inp, 'inp'],
    'vasp': [vasp.inp, ''],
    'espresso': [espresso.inp, 'inp'],
    'abinit': [abinit.inp, 'inp'],
    'nwchem': [nwchem.inp, 'inp'],
    'gaussian': [gaussian.inp, 'com'],
    'horton': [horton.inp, 'inp'],
    'ofdft': [ofdft.inp, 'inp'],
    'bigdft': [bigdft.inp, 'yaml'],
  }

  if type(molecule) is str:
    molecule = qtk.Molecule(molecule, **kwargs)

  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.qmcode
  p_str = kwargs['program'].lower()
  p = inp_dict[p_str][0]
  if 'charge' in kwargs:
    molecule.charge = kwargs['charge']
  kwargs['extension'] = inp_dict[p_str][1]
  return p(molecule, **kwargs)

def QMOut(out=None, **kwargs):
  out_dict = {
    'cpmd': cpmd.out,
    'vasp': vasp.out,
    'abinit': abinit.out,
    'espresso': espresso.out,
    'nwchem': nwchem.out,
    'horton': horton.out,
    'gaussian': gaussian.out,
    'bigdft': bigdft.out,
  }


  if 'unit' in kwargs:
    unit = kwargs['unit']
  else:
    unit = 'Ha'

  output = qtk.QMOutput()

  if out is not None and not os.path.exists(out):
    qtk.warning("%s not found" % out)
    #return qtk.QMOut()
  else:
    if 'program' in kwargs:
      p_str = kwargs['program']
      if 'debug' in kwargs and kwargs['debug']:
        #return out_dict[p_str](out).inUnit(unit)
        output = out_dict[p_str](out, **kwargs)
      else:
        try:
          #return out_dict[p_str](out)
          output = out_dict[p_str](out, **kwargs)
        except Exception as e:
          qtk.warning("%s failed with message: %s" % (out, e))
          #qout = qtk.QMOut()
          output.path, output.name = os.path.split(out)
          #return output.inUnit(unit)
    else:
      for p in out_dict.itervalues():
        try:
          output = p(out)
        except:
          pass
      qtk.warning("something wrong with output file, "+\
                  "pass 'program' eplicitly")
  return output.inUnit(unit)
