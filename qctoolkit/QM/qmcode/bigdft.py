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
    self.backup()

  def run(self, name=None, **kwargs):
    self.setting.update(kwargs)
    return univ.runCode(self, WaveletInput, name, **self.setting)

  def write(self, name=None, **kwargs):
    if 'no_reset' not in kwargs or not kwargs['no_reset']:
      self.reset()
    self.setting.update(kwargs)
    self.setting['yaml'] = True
    self.setting['root_dir'] = name
    inp, molecule = \
      super(WaveletInput, self).write(name, **self.setting)

    dft_list = {
                 'lda': 1,
                 'pbe': 11,
                 'pbexc': -101130,
                 'pbesol': -116133,
                 'pbe0': -406,
                 'blyp': -106131,
                 'b3lyp': -402,
               }
    
    if self.setting['theory'] in dft_list:
      theory = dft_list[self.setting['theory']]
    else:
      qtk.exit("theory: %s not implemented" % self.setting['theory'])
 
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

    self.pp_path = None
    if 'pp_path' not in self.setting:
      self.pp_path = qtk.setting.bigdft_pp
    else:
      self.pp_path = self.setting['pp_path']

    pp_base = ['pbe']
    if 'pp_theory' not in self.setting:
      self.setting['pp_theory'] = self.setting['theory']
    pp_dir = None
    for base in pp_base:
      if base in self.setting['pp_theory']:
        pp_dir = os.path.join(self.pp_path, base)
    if not pp_dir:
      qtk.warning('PP is not implemented for theory:%s. '+\
                  'To include, set pp_theory=implemented_theory,'+\
                  ' e.g. pbe')
      
    dft = {
            'hgrids': 'fast',
            'rmult': yList([3.5, 9.0]),
            'nrepmax': 'accurate',
            'disablesym': 'Yes',
            'ixc': theory,
          }

    posinp = {
               'units': 'angstroem',
             }
    positions = []
    pp_files = []
    posinp['positions'] = positions
    for i in range(molecule.N):
      entry = {molecule.type_list[i]: yList(list(molecule.R[i]))}
      positions.append(entry)
      pp_file = os.path.join(pp_dir, 
                             'psppar.' + molecule.type_list[i])
      pp_files.append(pp_file)

    if 'box' in self.setting and self.setting['box']:
      box = copy.deepcopy(self.setting['box'])
      if box[1] <= 0:
        box[1] = '.inf'
      posinp['cell'] = yList(box)
        

    data = {}
    data['dft'] = dft
    data['posinp'] = posinp

    content = reformat(yaml.dump(data, default_flow_style=False))

    inp.write(content)
    inp.close(dependent_files=pp_files)

    if name: return inp, name

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
    string = data[ind] + data[ind + 1]
    string = re.sub('\n', '', string)
    string = re.sub('.*{', '', string)
    string = re.sub('}.*', '', string)
    string = re.sub(' *', '', string)
    tmp = filter(None, string.split(','))
    self.scf_step = len(tmp)
    Et = {}
    for entry in tmp:
      name, value = entry.split(':')
      Et[name] = value
    self.detail = Et

    tmp = filter(lambda x: 'Energy (Hartree)' in x, data)
    self.Et = float(tmp[-1].split(':')[1])
