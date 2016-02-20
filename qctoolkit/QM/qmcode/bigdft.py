import qctoolkit as qtk
from qctoolkit.QM.wavelet_io import WaveletInput
from qctoolkit.QM.wavelet_io import WaveletOutput
import sys, os, re, copy, shutil
import qctoolkit.QM.qmjob as qmjob
import pkg_resources
import numpy as np
import universal as univ
import yaml

class inp(WaveletInput):
  """
  bigdft input class.
  """
  __doc__ = WaveletInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    WaveletInput.__init__(self, molecule, **kwargs)
    self.setting.update(**kwargs)

  def run(self, name=None, **kwargs):
    self.setting.update(kwargs)
    return univ.runCode(self, WaveletInput, name, **self.setting)

  def write(self, name=None, **kwargs):
    self.setting.update(kwargs)
    self.setting['yaml'] = True
    inp, molecule = \
      super(WaveletInput, self).write(name, **self.setting)
 
    def yList(data):
      _ = ' bracket '
      out = _
      for i in range(len(data)-1):
        out = out + str(data[i]) + ', '
      out = out + str(data[-1]) + _
      return out

    def reformat(content):
      out = re.sub("' bracket ", '[', content)
      out = re.sub(" bracket '", ']', out)
      out = re.sub("'", '', out)
      return out
      
    dft = {
            'hgrids': 'fast',
            'rmult': yList([3.5, 9.0]),
            'nrepmax': 'accurate',
            'output_wf': 1,
            'disablesym': 'Yes',
            'ixc': 11,
          }

    posinp = {
               'unit': 'angstroem',
             }
    positions = []
    posinp['position'] = positions
    for i in range(molecule.N):
      entry = {molecule.type_list[i]: yList(list(molecule.R[i]))}
      positions.append(entry)

    data = {}
    data['dft'] = dft
    data['posinp'] = posinp

    content = reformat(yaml.dump(data, default_flow_style=False))

    inp.write(content)
    inp.close()

    return univ.writeReturn(inp, name, **self.setting)

class out(WaveletOutput):
  def __init__(self, qmout, **kwargs):
    WaveletOutput.__init__(self, qmout, **kwargs)
    self.info = ''
    if qmout:
      self.getEt(qmout)

  def getEt(self, name):
    self.ino = os.path.splitext(name)[0]
    qmout = open(name)
    data = qmout.readlines()
    tmp = filter(lambda x: 'Energies' in x,  data)
    ind = data.index(tmp[-1])
    string = data[inp] + data[inp + 1]
    string = re.sub('\n', '', string)
    string = re.sub('.*{', '', string)
    string = re.sub('}.*', '', string)
    string = re.sub(' *', '', string)
    tmp = string.split(',')
    self.scf_step = len(tmp)
    Et = {}
    for entry in tmp:
      name, value = entry.split(':')
      Et[name] = value
    self.detail = Et

    tmp = filter(lambda x: 'Energy (Hartree)' in x, data)
    self.Et = float(tmp.split(':')[1])
