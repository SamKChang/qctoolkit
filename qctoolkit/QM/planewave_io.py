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
    angles = self.celldm[3:]
    if self.scale:
      lattice = [self.celldm[i]/self.scale[i] for i in range(3)]
    else:
      lattice = [self.celldm[i] for i in range(3)]
    lattice.extend(angles)
    fm = qtk.fractionalMatrix(lattice)
    self.setting['lattice'] = np.dot(fm, np.eye(3)).T

  def write(self, name=None, **kwargs):
    inp, molecule = \
      super(GenericQMInput, self).write(name, **setting)
    return inp, molecule

class PlanewaveOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
