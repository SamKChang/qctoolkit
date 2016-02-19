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

    data = {}
    data['name'] = 'value'
    dft = {}
    dft['ixc'] = 11
    dft['cutoff'] = 2
    data['dft'] = dft

    inp.write(yaml.dump(data, default_flow_style=False))
    inp.close()

    return univ.writeReturn(inp, name, **self.setting)

class out(WaveletOutput):
  def __init__(self, qmout, **kwargs):
    WaveletOutput.__init__(self, qmout, **kwargs)
    self.info = ''
    if qmout:
      self.getEt(qmout)

  def getEt(self, name):
    pass
