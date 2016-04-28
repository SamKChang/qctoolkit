import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
from veint import veint
from eeint import eeint
from neint import neint
from nnint import nnint
from vnint import vnint
from keint import keint
import numpy as np
from numpy import tensordot as td
import warnings 


def veMatrix(basis, coord, Z):
  basis_data, center, lm = basisData(basis)
  coord = np.array(coord) * 1.889725989
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return veint(basis_data, center, lm, coord, list(Z))

def vnMatrix(basis, coord, Z):
  basis_data, center, lm = basisData(basis, density=True)
  coord = np.array(coord) * 1.889725989
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return vnint(basis_data, center, lm, coord, list(Z))

def eeMatrix(basis):
  basis_data, center, lm = basisData(basis)
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return eeint(basis_data, center, lm)

def keMatrix(basis):
  basis_data, center, lm = basisData(basis)
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return keint(basis_data, center, lm)

def knMatrix(basis):
  basis_data, center, lm = basisData(basis)
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return knint(basis_data, center, lm)

def nnMatrix(basis):
  basis_data, center, lm = basisData(basis, density=True)
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return nnint(basis_data, center, lm)

def neMatrix(basis, fit_basis=None):
  if fit_basis is None:
    fit_basis = basis
  basis_data, center, lm = basisData(basis)
  fbasis_data, fcenter, flm = basisData(fit_basis, density=True)
  warnings.filterwarnings("ignore", category=DeprecationWarning) 
  return neint(basis_data, center, lm, fbasis_data, fcenter, flm)

def densityMatrix(qmout):
  """
  density matrix for closed shell system
  """
  psi_basis = qmout.basis
  occ = [i for i in range(qmout.n_ao)
         if qmout.occupation[i]==2][-1] + 1
  mo = qmout.mo_vectors
  mo_occ = mo[0:occ, :]
  return 2*np.dot(mo_occ.T, mo_occ)
  

def densityFitting(qmout, rho_basis=None):
  psi_basis = qmout.basis
  occ = [i for i in range(qmout.n_ao)
         if qmout.occupation[i]==2][-1] + 1
  D_matrix = densityMatrix(qmout)
  if rho_basis is None:
    rho_basis = psi_basis

  psi_basis_data, psi_center, psi_lm = basisData(psi_basis)
  rho_basis_data, rho_center, rho_lm = \
    basisData(rho_basis, density=True)
  NN = nnMatrix(rho_basis)
  NE = neMatrix(psi_basis, rho_basis)
  G = td(D_matrix, NE, axes=([0,1], [1,2]))
  
  return np.linalg.solve(NN, G)

def basisData(basis, **kwargs):
  if 'density' not in kwargs:
    kwargs['density'] = False
  centers = []
  exponents = []
  coefficients = []  
  lm_xyz = []
  n_gaussians = []
  for ao in basis:
    ng = 0
    for g in range(len(ao['exponents'])):
      if not kwargs['density']:
        exponents.append(ao['exponents'][g])
      else:
        exponents.append(2*ao['exponents'][g])
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
  if kwargs['density']:
    lm_xyz = 2*np.array(lm_xyz)
  else:
    lm_xyz = np.array(lm_xyz)
  return out, np.array(centers), lm_xyz
  

class GaussianBasisInput(GenericQMInput):
  """
  From GaussianBasis Input:
  generic class holder for gaussian basis qmcode. It provide basic
  default settings.
  ===
  """
  __doc__ = GenericQMInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if 'basis_set' not in kwargs:
      self.setting['basis_set'] = '6-31g'

class GaussianBasisOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
