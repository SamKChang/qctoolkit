import qctoolkit as qtk
from inp import GenericInput
import numpy as np

class AtomicBasisInput(GenericInput):
  def __init__(self, molecule, **kwargs):
    GenericInput.__init__(self, molecule, **kwargs)

