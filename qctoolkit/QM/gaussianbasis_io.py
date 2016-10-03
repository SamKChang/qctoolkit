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
import pkgutil
import os

from ofdft.libxc_dict import xc_dict
import ofdft.libxc_interface as xcio
ps_eggs_loader = pkgutil.find_loader('pyscf')
ps_found = ps_eggs_loader is not None
ht_eggs_loader = pkgutil.find_loader('horton')
ht_found = ht_eggs_loader is not None

selfPath = os.path.realpath(__file__)
selfPath = os.path.split(selfPath)[0] + 'ofdft'
xcpath = os.path.join(selfPath, 'libxc_exc.so')
xc_found = os.path.exists(xcpath)
if ps_found:
  from pyscf import gto
  import pyscf as ps
else:
  pass
if ht_found:
  from horton import BeckeMolGrid
else:
  pass

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
      self.setting['basis_set'] = 'def2-tzvp'

class GaussianBasisOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)

  def mo_g09_nwchem(self):
    mo = self.mo_vectors
    ind = np.arange(len(mo[0]))
    itr = 0

    order = {
     'd': [0, 3, 4, 1, 5, 2],
     'f': [0, 4, 5, 3, 9, 6, 1, 8, 7, 2],
    }

    while itr < len(mo[0]):
      bStr = self.basis[itr]['type']
      if len(bStr) > 2:
        key = bStr[0]
        if key in order.keys():
          ind_itr = ind[itr: itr+len(order[key])]
          ind_itr = ind_itr[order[key]]
          for i in range(len(order[bStr[0]])):
            ind[itr + i] = ind_itr[i]
          itr += len(order[bStr[0]]) - 1
        else:
          qtk.exit("basis reordering for %s orbital " % key\
                 + "not implemented yet")
      itr += 1
    return mo[:, ind]

  def basisFormat():
    if not hasattr(self, 'basis'):
      qtk.exit("no basis found")
    pass

  def getBeckeGrid(self, grid='fine'):
    """
    coarse, medium, fine, veryfine, ultrafine and insane
    """
    if not hasattr(self, 'molecule'):
      qtk.exit("no molecule structure found")
    molecule = self.molecule
    coord = np.array(np.atleast_2d(molecule.R*1.8897261245650618))
    self.grid = BeckeMolGrid(coord, 
                             molecule.Z.astype(int), 
                             molecule.Z,
                             grid)

    mol_str = []
    for i in range(molecule.N):
      atm_str = [molecule.type_list[i]]
      for j in range(3):
         atm_str.append(str(molecule.R[i,j]))
      mol_str.append(' '.join(atm_str))
    mol_str = '; '.join(mol_str)

    mol = gto.Mole()
    #mol.build(atom=mol_str, basis=self.setting['basis_set'])
    if hasattr(self, 'basis_name'):
      basis = self.basis_name
    mol.build(atom=mol_str, basis=basis)
    self.mol = mol

  def getPhi(self, cartesian=True, grid='fine'):
    self.getBeckeGrid(grid)
    coords = self.grid.points
    if cartesian:
      mode = "GTOval_cart"
    else:
      mode = "GTOval_sph"
    self._phi = self.mol.eval_gto(mode, coords).T
    norm = np.dot(self._phi * self.grid.weights, self._phi.T)
    self._phi = self._phi / np.sqrt(np.diag(norm))[:, np.newaxis]
    return self._phi

  def getPsi(self, cartesian=True, grid='fine'):
    self.getPhi(cartesian, grid)
    if not hasattr(self, 'mo_vectors'):
      qtk.exit('mo_vectors not found')
    mo = self.mo_vectors
    if hasattr(self, 'program'):
      if self.program == 'gaussian':
        mo = self.mo_g09_nwchem()
    self._psi = np.dot(mo, self._phi)
    return self._psi

  def getRho(self, cartesian=True, grid='fine'):
    self.getPsi(cartesian, grid)
    if not hasattr(self, 'occupation'):
      qtk.exit("occupation number not found")
    self.rho = np.zeros(self.grid.size)
    for i in range(len(self.occupation)):
      n = self.occupation[i]
      psi = self._psi[i]
      self.rho += n * psi**2
    return self.rho

  def getDipole(self, cartesian=True, grid='find', unit='debye'):
    if not hasattr(self, 'molecule'):
      qtk.exit('molecule structure not found')
    if not hasattr(self, 'rho'):
      self.getRho(cartesian, grid)
    pQ = np.array(
      [sum(self.molecule.Z * self.molecule.R[:,i]) for i in range(3)]
    )
    pQ = pQ * 1.8897259885789
    pq = np.array(
      [self.grid.integrate(self.rho * self.grid.points[:,i]) 
       for i in range(3)]
    )
    pq = pq 
    mu = pQ - pq
    if unit == 'debye':
      return mu / 0.393430307
    else:
      return mu

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
  

