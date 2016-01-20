import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
from veint import veint
from eeint import eeint
import numpy as np

def veMatrix(basis, coord, Z):
  basis_data, center, lm = basisData(basis)
  coord = np.array(coord) * 1.889725989
  return veint(basis_data, center, lm, coord, list(Z))

def eeMatrix(basis):
  basis_data, center, lm = basisData(basis)
  return eeint(basis_data, center, lm)

def keMatrix(basis):
  pass

def nnMatrix(basis):
  pass

def neMatrix(basis):

def basisData(basis):
  centers = []
  exponents = []
  coefficients = []  
  lm_xyz = []
  n_gaussians = []
  for ao in basis:
    ng = 0
    for g in range(len(ao['exponents'])):
      exponents.append(ao['exponents'][g])
      coefficients.append(ao['coefficients'][g])
      ng = ng + 1
      if g == len(ao['exponents'])-1:
        n_gaussians.append(g+1)
        centers.append(list(ao['center']))
        lm_xyz.append([ao['type'].count('x'), 
                       ao['type'].count('y'), 
                       ao['type'].count('z')])
  out = {}
  out['exponents'] = exponents
  out['coefficients'] = coefficients
  out['n_gaussians'] = n_gaussians
  return out, np.array(centers), np.array(lm_xyz)
  

class AtomicBasisInput(GenericQMInput):
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if 'basis_set' not in kwargs:
      self.setting['basis_set'] = '6-31g'

class AtomicBasisOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
