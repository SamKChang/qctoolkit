import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
from veint import veint
from eeint import eeint
from eekernel import eekernel
from neint import neint
from nnint import nnint
from vnint import vnint
from keint import keint
import numpy as np
from numpy import tensordot as td
import warnings 
import pkgutil
import os
import copy

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
  try:
    from horton import BeckeMolGrid
  except:
    pass
else:
  pass
if xc_found:
  try:
    from ofdft.libxc_exc import libxc_exc
    from ofdft.libxc_dict import xc_dict
  except:
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

  def copy(self):
    if hasattr(self, 'mol'):
      del self.mol
    if hasattr(self, 'grid'):
      del self.grid
    return copy.deepcopy(self)

  def mo_g09_nwchem(self):
    mo = self.mo_vectors
    ind = np.arange(len(mo[0]))
    itr = 0

    # hard coded reordering for d and f orbitals
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

  def basisFormat(self):

    l_dict = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5}

    if not hasattr(self, 'basis'):
      qtk.exit("no basis found")
    basis = {}
    for b in self.basis:
      if b['atom'] not in basis:
        basis[b['atom']] = []
      l = None
      for k, v in l_dict.items():
        if k in b['type']:
          l = v
      if l is None:
        qtk.exit('basis extraction failed at' + str(b))
      cef = b['coefficients']
      exp = b['exponents']
      ncef = np.asarray(cef)
      nexp = np.asarray(exp)
      contained = False
      for data in basis[b['atom']]:
        if data[0] == l:
          data_c = []
          data_e = []
          for d in data[1:]:
            data_c.append(d[1])
            data_e.append(d[0])
          try:
            norm_c = np.linalg.norm(ncef - data_c)
            norm_e = np.linalg.norm(nexp - data_e)
            if norm_e < 1E-5 and norm_c < 1E-5:
              contained = True
              break
          except:
            pass
      if not contained:
        new_b = [list(ec) for ec in zip(exp, cef)]
        b_data = [l]
        b_data.extend(new_b)
        basis[b['atom']].append(b_data)
    self.pybasis = basis

  def getBeckeGrid(self, resolution='fine', new=False, **kwargs):
    """
    coarse, medium, fine, veryfine, ultrafine and insane
    """
    if not hasattr(self, 'molecule'):
      qtk.exit("no molecule structure found")
    if new or not hasattr(self, 'grid'):
      molecule = self.molecule
      coord = np.array(np.atleast_2d(molecule.R*1.8897261245650618))
      self.grid = BeckeMolGrid(coord, 
                               molecule.Z.astype(int), 
                               molecule.Z,
                               resolution)
  
      mol_str = []
      for i in range(molecule.N):
        atm_str = [molecule.type_list[i]]
        for j in range(3):
           atm_str.append(str(molecule.R[i,j]))
        mol_str.append(' '.join(atm_str))
      mol_str = '; '.join(mol_str)
  
      if 'gto_kwargs' in kwargs:
        mol = gto.Mole(**kwargs['gto_kwargs'])
      else:
        mol = gto.Mole()

      #mol.build(atom=mol_str, basis=self.setting['basis_set'])
      if hasattr(self, 'basis_name'):
        basis = self.basis_name
      self.basisFormat()
      mol.build(atom=mol_str, basis=self.pybasis)
      self.mol = mol
      del_list = ['_phi', '_psi', '_dphi', '_dpsi', '_rho', '_drho']
      for p in del_list:
        if hasattr(self, p):
          delattr(self, p)

  def getPhi(self, cartesian=True, resolution='fine', new=False, **kwargs):
    if 'gridpoints' in kwargs:
      new = True
    if new or not hasattr(self, '_phi'):
      self.getBeckeGrid(resolution, new, **kwargs)
      coords = self.grid.points
      if 'gridpoints' not in kwargs:
        grid_coords = None
      else:
        grid_coords = np.array(kwargs['gridpoints']).astype(float)
      if cartesian:
        mode = "GTOval_cart"
      else:
        mode = "GTOval_sph"
      self._phi = self.mol.eval_gto(mode, coords).T
      norm = np.dot(self._phi * self.grid.weights, self._phi.T)

      if grid_coords is not None:
        self._phi = self.mol.eval_gto(mode, grid_coords).T
        self._phi = self._phi / np.sqrt(np.diag(norm))[:, np.newaxis]
      else:
        self._phi = self._phi / np.sqrt(np.diag(norm))[:, np.newaxis]
    return self._phi

  def getDPhi(self, cartesian=True, resolution='fine', new=False, **kwargs):
    if 'gridpoints' in kwargs:
      new = True
    if new or not hasattr(self, '_dphi'):
      self.getBeckeGrid(resolution, new, **kwargs)
      if 'gridpoints' not in kwargs:
        coords = self.grid.points
      else:
        coords = np.array(kwargs['gridpoints']).astype(float)
      if cartesian:
        mode = "GTOval_ip_cart"
      else:
        mode = "GTOval_ip_sph"
      self._dphi = self.mol.eval_gto(mode, coords, comp=3).T
      
    return self._dphi

  def getPsi(self, cartesian=True, resolution='fine', new=False, **kwargs):
    if 'gridpoints' in kwargs:
      new = True
    if new or not hasattr(self, '_psi'):
      self.getPhi(cartesian, resolution, new, **kwargs)
      if not hasattr(self, 'mo_vectors'):
        qtk.exit('mo_vectors not found')
      mo = self.mo_vectors
      if hasattr(self, 'program'):
        if self.program == 'gaussian':
          mo = self.mo_g09_nwchem()
      self._psi = np.dot(mo, self._phi)
    return self._psi

  def getDPsi(self, cartesian=True, resolution='fine', new=False, **kwargs):
    if 'gridpoints' in kwargs:
      new = True
    if new or not hasattr(self, '_dpsi'):
      self.getDPhi(cartesian, resolution, new, **kwargs)
      if not hasattr(self, 'mo_vectors'):
        qtk.exit('mo_vectors not found')
      mo = self.mo_vectors
      if hasattr(self, 'program'):
        if self.program == 'gaussian':
          mo = self.mo_g09_nwchem()
      self._dpsi = np.dot(mo, np.swapaxes(self._dphi, 0, 1))
    return self._dpsi

  def getRho(self, cartesian=True, resolution='fine', new=False, occupation=None, **kwargs):
    if 'gridpoints' in kwargs:
      new = True
    if new or not hasattr(self, '_rho'):
      self.getPsi(cartesian, resolution, new, **kwargs)
      if not hasattr(self, 'occupation'):
        qtk.exit("occupation number not found")
      occ = np.array(self.occupation)
      if occupation is not None:
        for i in range(len(occupation)):
          occ[i] = occupation[i]
      self._rho = np.sum(self._psi**2 * occ[:, np.newaxis], axis = 0)
    return self._rho

  def getDRho(self, cartesian=True, resolution='fine', new=False, occupation=None, **kwargs):
    if 'gridpoints' in kwargs:
      new = True
    if new or not hasattr(self, '_drho'):
      if not hasattr(self, '_psi'):
        _psi = self.getPsi(cartesian, resolution, new, **kwargs)
      if not hasattr(self, '_dpsi'):
        _dpsi = self.getDPsi(cartesian, resolution, new, **kwargs)
      if not hasattr(self, 'occupation'):
        qtk.exit("occupation number not found")
      self._psi, self._dpsi = _psi, _dpsi
      if occupation is None:
        occ = np.array(self.occupation)
      else:
        occ = np.array(occupation)
      self._drho = 2 * np.sum(
        self._psi[..., np.newaxis] * self._dpsi \
        * occ[:, np.newaxis, np.newaxis],
        axis = 0
      )
    return self._drho

  def _freeSetup(self, coord, gridpoints, **kwargs):
    assert self.molecule.N == 1
    new = self.copy()
    new.molecule.R[0] = np.array(coord) / 1.8897261245650618
    if 'spin' not in kwargs:
      if new.molecule.getValenceElectrons() % 2 != 0:
        kw = {'gto_kwargs': {'spin': 1}}
      else:
        kw = {}
    else:
      kw = {'gto_kwargs': {'spin': kwargs['spin']}}
    if type(gridpoints[0][0]) is not float:
      pl = np.array(gridpoints).astype(float)
    else:
      pl = gridpoints
    return new, kw, pl

  def freePsi(self, coord, gridpoints, **kwargs):
    new, kw, pl = self._freeSetup(coord, gridpoints)
    return new.getRho(gridpoints = pl, **kw)

  def freeDPsi(self, coord, gridpoints, **kwargs):
    new, kw, pl = self._freeSetup(coord, gridpoints)
    return new.getDRho(gridpoints = pl, **kw)

  def freeRho(self, coord, gridpoints, **kwargs):
    new, kw, pl = self._freeSetup(coord, gridpoints)
    return new.getRho(gridpoints = pl, **kw)

  def freeDRho(self, coord, gridpoints, **kwargs):
    new, kw, pl = self._freeSetup(coord, gridpoints)
    return new.getDRho(gridpoints = pl, **kw)

  def _e_setting(self, gridpoints, **kwargs):
    new = self.copy()
    kw = {}
    if 'spin' not in kwargs:
      if new.molecule.getValenceElectrons() % 2 != 0:
        kw['gto_kwargs'] = {'spin': 1}
    else:
      kw['gto_kwargs'] = {'spin': spin}
    if 'resolution' in kwargs:
      kw['resolution'] = kwargs['resolution']
    if 'cartesian' in kwargs:
      kw['cartesian'] = kwargs['cartesian']
    if 'free_coord' in kwargs:
      assert self.molecule.N == 1
      # assume bohr input
      coord = kwargs['free_coord']
      new.molecule.R[0] = np.array(coord) / 1.8897261245650618
      for b in new.basis:
        b['center'] = np.array(coord)
    occ = np.array(new.occupation)
    gridpoints = np.atleast_2d(gridpoints).astype(float)
    return new, occ, kw, gridpoints

  def e_kin(self, gridpoints, **kwargs):
    """
    kinetic energy density
    """
    new, occ, kw, gridpoints = self._e_setting(gridpoints, **kwargs)
    dpsi = new.getDPsi(gridpoints = gridpoints, **kw)
    dpsi2 = np.linalg.norm(dpsi, axis=-1) ** 2
    return (0.5 * dpsi2 * occ[:, np.newaxis]).sum(axis=0)

  def e_ext(self, gridpoints, **kwargs):
    """
    external potential energy density
    """
    new, occ, kw, gridpoints = self._e_setting(gridpoints, **kwargs)

    ext = np.zeros(len(gridpoints))
    for I in range(new.molecule.N):
      RI = new.molecule.R[I] * 1.8897261245650618
      ZI = new.molecule.Z[I]
      R_list = gridpoints - RI
      ext -= ZI / np.linalg.norm(R_list, axis=1)

    rho = new.getRho(gridpoints = gridpoints, **kw)
    #return (psi * occ[:, np.newaxis]).sum(axis=0) * ext
    return rho * ext

  def e_xc(self, gridpoints, dft='pbe', **kwargs):
    new, occ, kw, gridpoints = self._e_setting(gridpoints, **kwargs)
    xc_list = [xc_dict[f] for f in qtk.setting.libxc_dict[dft]]

    drho = new.getDRho(gridpoints=gridpoints, new=True, **kw)
    sigma = np.sum(drho ** 2, axis = 1)
    rho = new.getRho(gridpoints = gridpoints, **kw)
    x = libxc_exc(rho, sigma, len(rho), xc_list[0])
    c = libxc_exc(rho, sigma, len(rho), xc_list[1])
    return x, c

  def e_coulomb(self, gridpoints, **kwargs):
    new, occ, kw, gridpoints = self._e_setting(gridpoints, **kwargs)

    if 'batch_size' in kwargs:
      batch_size = kwargs['batch_size']
      del kwargs['batch_size']
    else:
      batch_size = 1

    out = []
    itr = 1
    coord = gridpoints
    for chunk in np.array_split(coord, batch_size):
      qtk.progress("e_coulomb", "processing batch: %d\n" % itr)
      itr += 1
      kernel = new.coulombKernel(
        chunk / 1.8897261245650618, 
        **kwargs
      )
      rho = new.getRho(gridpoints = chunk, **kw)
      out.append(0.5 * rho * kernel)
    return np.concatenate(out)

  def epsilon(self, gridpoints, dft='pbe', **kwargs):
    new, occ, kw, gridpoints = self._e_setting(gridpoints, **kwargs)
    kin = self.e_kin(gridpoints, **kw)
    ext = self.e_ext(gridpoints, **kw)
    x, c = self.e_xc(gridpoints, dft=dft, **kw)
    coulomb = self.e_coulomb(gridpoints, **kwargs)
    return kin + ext + x + c + coulomb

  def getDipole(self, cartesian=True, resolution='fine', unit='debye', 
                component='full'):
    if not hasattr(self, 'molecule'):
      qtk.exit('molecule structure not found')
    if not hasattr(self, '_rho'):
      self.getRho(cartesian, resolution)
    pQ = np.array(
      [sum(self.molecule.Z * self.molecule.R[:,i]) for i in range(3)]
    )
    pQ = pQ * 1.8897259885789

    if component in ['full', 'ele']:
      pq = np.array(
        [self.grid.integrate(self._rho * self.grid.points[:,i]) 
         for i in range(3)]
      )

    if component == 'full':
      mu = pQ - pq
    elif component == 'ele':
      mu = pq
    elif component == 'nuc':
      mu = pQ
    if unit == 'debye':
      return mu / 0.393430307
    else:
      return mu

  def keMatrix(self):
    return keMatrix(self.basis)

  def densityMatrix(self):
    return densityMatrix(self)

  def veMatrix(self, coord = None, Z = None):
    if coord is None:
      coord = self.molecule.R
    if Z is None:
      Z = self.molecule.Z
    return veMatrix(self.basis, coord, Z)

  def eeKernel(self, coord=None, batch_size=1, **kwargs):
    if coord is None:
      coord = self.molecule.R
    else:
      coord = np.atleast_2d(coord).astype(float)
    out = []
    itr = 1
    for chunk in np.array_split(coord, batch_size):
      if itr > 1:
        qtk.progress("eeKernel", "processing batch: %d\n" % itr)
      itr += 1
      out.append(eeKernel(self.basis, chunk))
    return np.concatenate(out)

  def coulombKernel(self, coord = None, **kwargs):
    k = self.eeKernel(coord, **kwargs)
    mo = self.mo_vectors
    # out = np.diagonal(
    #   np.einsum('is,jl,kls->kij', mo, mo, k),
    #   axis1=1, axis2=2,
    # ).sum(1)
    # 
    # diagonal: np.einsum('is,il, kls->ki', mo, mo, k)
    #   or np.einsum('ii->i', a) for diagonal
    # sum diagonal: np.einsum('is,il, kls->k', mo, mo, k)
    #   or np.einsum('ii', a) for trace
    out = np.einsum('is,il, kls->k', mo, mo, k)
    return out

  def eeMatrix(self):
    if not hasattr(self, '_eeMatrix'):
      self._eeMatrix = eeMatrix(self.basis)
    return self._eeMatrix

  def EJ(self):
    ee = 2*self.eeMatrix()
    td = np.tensordot
    mo = self.mo_vectors
    out = td(mo, ee, axes=(1, 0))
    out = td(mo, out, axes=(1, 1))
    out = td(mo, out, axes=(1, 2))
    out = td(mo, out, axes=(1, 3))
    occ = int(np.sum(self.occupation) / 2.)
    EJ = [out[a,a,b,b] for a in range(occ) for b in range(occ)]
    return sum(EJ)

  def EX(self):
    ex = -np.swapaxes(self.eeMatrix(), 1,2)
    td = np.tensordot
    mo = self.mo_vectors
    out = td(mo, ex, axes=(1, 0))
    out = td(mo, out, axes=(1, 1))
    out = td(mo, out, axes=(1, 2))
    out = td(mo, out, axes=(1, 3))
    occ = int(np.sum(self.occupation) / 2.)
    EX = [out[a,a,b,b] for a in range(occ) for b in range(occ)]
    return sum(EX)

  def EK(self):
    km = self.keMatrix()
    dm = self.densityMatrix()
    return np.trace(np.dot(dm, km))

  def Eext(self):
    vext = self.veMatrix()
    dm = self.densityMatrix()
    return np.trace(np.dot(dm, vext))

  # density fitting results not reliable... yet!
  #def vnMatrix(self, coord = None, Z = None):
  #  if coord is None:
  #    coord = self.molecule.R
  #  if Z is None:
  #    Z = self.molecule.Z
  #  return vnMatrix(self.basis, coord, Z)

  #def knMatrix(self):
  #  return knMatrix(self.basis)

  #def nnMatrix(self):
  #  return nnMatrix(self.basis)

  #def neMatrix(self):
  #  return neMatrix(self.basis)

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

def eeKernel(basis, coord):
  basis_data, center, lm = basisData(basis)
  coord = np.array(coord) * 1.889725989
  warnings.filterwarnings("ignore", category=DeprecationWarning)
  dummy_Z = [1. for _ in range(len(coord))]
  return eekernel(basis_data, center, lm, coord, dummy_Z).swapaxes(0,-1)

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
  

