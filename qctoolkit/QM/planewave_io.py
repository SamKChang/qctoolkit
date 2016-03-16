import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
import universal as univ
import numpy as np

class PlanewaveInput(GenericQMInput):
  """
  From PlanwaveInput:
  generic class holder for plane wave qmcode. It provide basic
  default settings.
  """
  __doc__ = GenericQMInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if 'cutoff' not in kwargs:
      self.setting['cutoff'] = 100
    if not self.setting['periodic'] and 'isolation' not in kwargs:
      self.setting['isolation'] = 'mt'
    self.pp_files = []
    if 'periodic' in self.setting and self.setting['periodic']:
      self.celldm2lattice()

    univ.getCelldm(self) 

  def celldm2lattice(self):
    cd = self.setting['celldm']
    if 'scale' in self.setting:
      sc = self.setting['scale']
    else:
      sc = [1.0 for i in range(3)]
    self.setting['lattice'] = qtk.celldm2lattice(cd, scale=sc)

  def write(self, name=None, **kwargs):
    if self.setting['periodic']:
      self.celldm2lattice()
    inp, molecule = \
      super(GenericQMInput, self).write(name, **setting)
    return inp, molecule

class PlanewaveOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
