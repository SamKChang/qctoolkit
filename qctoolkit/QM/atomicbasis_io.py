import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
from veint import veint
from eeint import eeint
from neint import neint
from nnint import nnint
from vnint import vnint
import numpy as np
from numpy import tensordot as td

def veMatrix(basis, coord, Z):
  basis_data, center, lm = basisData(basis)
  coord = np.array(coord) * 1.889725989
  return veint(basis_data, center, lm, coord, list(Z))

def vnMatrix(basis, coord, Z):
  basis_data, center, lm = basisData(basis)
  coord = np.array(coord) * 1.889725989
  return vnint(basis_data, center, lm, coord, list(Z))

def eeMatrix(basis):
  basis_data, center, lm = basisData(basis)
  return eeint(basis_data, center, lm)

def keMatrix(basis):
  pass

def nnMatrix(basis):
  basis_data, center, lm = basisData(basis)
  return nnint(basis_data, center, lm)

def neMatrix(basis, fit_basis=None):
  if fit_basis is None:
    fit_basis = basis
  basis_data, center, lm = basisData(basis)
  fbasis_data, fcenter, flm = basisData(fit_basis)
  return neint(basis_data, center, lm, fbasis_data, fcenter, flm)

def densityMatrix(qmout):
  psi_basis = qmout.basis
  occ = [i for i in range(qmout.n_ao)
         if qmout.occupation[i]==2][-1] + 1
  mo = qmout.mo_vectors
  mo_occ = mo[0:occ, :]
  return np.outer(mo_occ, mo_occ)
  

def densityFitting(qmout, rho_basis=None):
  psi_basis = qmout.basis
  occ = [i for i in range(qmout.n_ao)
         if qmout.occupation[i]==2][-1] + 1
  D_matrix = densityMatrix(qmout)
  if rho_basis is None:
    rho_basis = psi_basis

  psi_basis_data, psi_center, psi_lm = basisData(psi_basis)
  rho_basis_data, rho_center, rho_lm = basisData(rho_basis)
  NN = nnMatrix(rho_basis)
  NE = neMatrix(psi_basis, rho_basis)
  G = td(D_matrix, NE, axes=([0,1], [1,2]))
  
  return np.linalg.solve(NN, G)

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
