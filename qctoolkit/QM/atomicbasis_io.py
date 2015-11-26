import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
from gcint import gcint
import numpy as np

def gcMatrix(basis, coord, Z):
  basis_data, center, lm = basisData(basis)
  return gcint(basis_data, center, lm, coord, list(Z))

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
        if ao['type'] == 's':
          lm_xyz.append([0,0,0])
        elif ao['type'] == 'px':
          lm_xyz.append([1,0,0])
        elif ao['type'] == 'py':
          lm_xyz.append([0,1,0])
        elif ao['type'] == 'pz':
          lm_xyz.append([0,0,1])
        elif ao['type'] == 'dxx':
          lm_xyz.append([2,0,0])
        elif ao['type'] == 'dxy':
          lm_xyz.append([1,1,0])
        elif ao['type'] == 'dxz':
          lm_xyz.append([1,0,0])
        elif ao['type'] == 'dyy':
          lm_xyz.append([0,2,0])
        elif ao['type'] == 'dyz':
          lm_xyz.append([0,1,1])
        elif ao['type'] == 'dzz':
          lm_xyz.append([0,0,2])
        elif ao['type'] == 'fxxx':
          lm_xyz.append([3,0,0])
        elif ao['type'] == 'fxxy':
          lm_xyz.append([2,1,0])
        elif ao['type'] == 'fxxz':
          lm_xyz.append([2,0,1])
        elif ao['type'] == 'fxyy':
          lm_xyz.append([1,2,0])
        elif ao['type'] == 'fxyz':
          lm_xyz.append([1,1,1])
        elif ao['type'] == 'fxzz':
          lm_xyz.append([1,0,2])
        elif ao['type'] == 'fyyy':
          lm_xyz.append([0,3,0])
        elif ao['type'] == 'fyyz':
          lm_xyz.append([0,2,1])
        elif ao['type'] == 'fyzz':
          lm_xyz.append([0,1,2])
        elif ao['type'] == 'fzzz':
          lm_xyz.append([0,0,3])
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
    GenericQMOutput.__init__(self, output=None, **kwargs)
