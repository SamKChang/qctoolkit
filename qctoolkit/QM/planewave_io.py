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

    if 'periodic' not in kwargs:
      self.setting['periodic'] = True
    if 'symmetry' not in kwargs and self.setting['periodic']:
      self.setting['symmetry'] = 'orthorhombic'
    if 'cutoff' not in kwargs:
      self.setting['cutoff'] = 100
    if not self.setting['periodic'] and 'isolation' not in kwargs:
      self.setting['isolation'] = 'mt'
    self.pp_files = []
    self.pp_path = ''

    univ.getCelldm(self) 

  def celldm2lattice(self):
    cd = self.setting['celldm']
    if cd[3]==0 and cd[4]==0 and cd[5]==0:
      self.setting['lattice'] = np.array([[cd[0],   0.0,   0.0],
                                          [  0.0, cd[1],   0.0],
                                          [  0.0,   0.0, cd[2]]])
  def center(self, coord):
    self.molecule.center(coord)
  def shift(self, coord):
    self.molecule.shift(coord)

  def write(self, name=None, **kwargs):
    inp, molecule = \
      super(GenericQMInput, self).write(name, **setting)
    return inp, molecule

class PlanewaveOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
