import qctoolkit as qtk
from inp import GenericQMInput
import numpy as np

class AtomicBasisInput(GenericQMInput):
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

