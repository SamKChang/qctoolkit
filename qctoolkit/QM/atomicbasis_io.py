import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
import numpy as np

class AtomicBasisInput(GenericQMInput):
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if 'basis_set' not in kwargs:
      self.setting['basis_set'] = '6-31g'

class AtomicBasisOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output=None, **kwargs)
