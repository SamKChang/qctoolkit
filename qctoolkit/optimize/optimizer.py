from qctoolkit import *
from ccs import *

class CCSOptimizer(object):
  def __init__(self, target_function, structure, ccs_parameters):
    self._ccs_space = MoleculeSpan(structure, ccs_parameters)
