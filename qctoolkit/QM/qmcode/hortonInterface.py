import qctoolkit as qtk
import pkgutil
eggs_loader = pkgutil.find_loader('horton')
found = eggs_loader is not None
if found:
  try:
    from horton import *
  except:
    pass
else:
  pass
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import periodictable as pt

class inp(GaussianBasisInput):
  """
  horton input class. 
  """
  __doc__ = GaussianBasisInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):

    if 'theory' not in kwargs:
      kwargs['theory'] = 'hf'

    if 'save_c_type' not in kwargs:
      kwargs['save_c_type'] = True

    if not found:
      qtk.exit("horton module not found.")
    if 'wf_convergence' not in kwargs:
      kwargs['wf_convergence'] = 1e-06
    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()

    coord = molecule.R * qtk.setting.a2b

    mol = IOData(coordinates=coord, numbers=molecule.Z)
    obasis = get_gobasis(mol.coordinates, mol.numbers,
                         self.setting['basis_set'])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, 
                        mol.pseudo_numbers)

    olp = obasis.compute_overlap()
    kin = obasis.compute_kinetic()
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    er = obasis.compute_electron_repulsion()

    exp_alpha = orb_alpha = Orbitals(obasis.nbasis)
    
    # Initial guess
    guess_core_hamiltonian(olp, kin + na, exp_alpha)

    external = {'nn': compute_nucnuc(mol.coordinates, 
                                     mol.pseudo_numbers)}
    if self.setting['theory'] == 'hf':
      terms = [
          RTwoIndexTerm(kin, 'kin'),
          RDirectTerm(er, 'hartree'),
          RExchangeTerm(er, 'x_hf'),
          RTwoIndexTerm(na, 'ne'),
      ]
    else:
      #libxc_term = RLibXCHybridGGA('xc_%s' % self.setting['theory'])
      libxc_term = RLibXCHybridGGA('xc_pbeh')
      terms = [
          RTwoIndexTerm(kin, 'kin'),
          RGridGroup(obasis, grid, [libxc_term]),
          RExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
          RExchangeTerm(er, 'x_hf'),
          RTwoIndexTerm(na, 'ne'),
      ]
    ham = REffHam(terms, external)

    occ_model = AufbauOccModel(
      int((sum(self.molecule.Z) - self.molecule.charge)/ 2.)
    )

    occ = np.zeros(olp.shape[0])
    N = int(sum(self.molecule.Z) - self.molecule.charge)
    occ[:N/2] = 1
    C = copy.deepcopy(exp_alpha.coeffs.__array__())
    dm = (C * occ).dot(C.T)

    if self.setting['save_c_type']:
      self.ht_mol = mol
      self.ht_grid = grid
      self.ht_external = external
      self.ht_obasis = obasis
      self.ht_occ_model = occ_model
      self.ht_olp = olp
      self.ht_kin = kin
      self.ht_na = na
      self.ht_er = er
      self.ht_exp_alpha = exp_alpha
      self.ht_terms = terms
      #self.ht_dm_alpha = dm_alpha
      self.ht_ham = ham
    else:
      self.grid_points = grid.points
      self.grid_weights = grid.weights

    self.olp = olp.__array__()
    self.kin = kin.__array__()
    self.nn = external['nn']
    try:
      d, U = np.linalg.eigh(self.olp)
      self.X = U.dot(np.diag(np.sqrt(1/d)).dot(U.T))
    except Exception as err:
      qtk.warning(err)
    self.na = na.__array__()
    self.er = er.__array__()
    self.mov = exp_alpha.coeffs.__array__()
    self.initial_mov = C
    self.initial_dm = dm
    self.occ = occ
    self.occupation = 2*occ
    self.v_ext = self.na
    #self.dm = dm_alpha.__array__()

  def run(self, name=None, **kwargs):
    scf_solver = PlainSCFSolver(threshold=self.setting['wf_convergence'], maxiter=self.setting['scf_step'])
    #scf_solver = CDIISSCFSolver(1e-6)
    scf_solver(
      self.ht_ham, 
      self.ht_olp, 
      self.ht_occ_model, 
      self.ht_exp_alpha
    )

    self.ht_solver = scf_solver

    try:
      self.Et = self.ht_ham.cache['energy']
      self.mov = self.ht_exp_alpha.coeffs.__array__()
      self.mo_vectors = self.mov
      self.mo_eigenvalues = self.ht_exp_alpha.energies

    except Exception as err:
      qtk.warning('SCF did not converged: %s' % err)
      self.Et = np.nan
    self.mol = self.molecule

    out = qtk.QM.gaussianbasis_io.GaussianBasisOutput()

    for attr in dir(out):
      if not hasattr(self, attr):
        setattr(self, attr, getattr(out, attr))

    for attr in dir(self):
      if not attr.startswith('_'):
        setattr(out, attr, getattr(self, attr))

    return out
    
  def write(self, name=None, **kwargs):
    pass

  def dm(self, C=None):
    if C is None: C = self.mov
    occ = self.occ
    return (C * occ).dot(C.T)

  def getRho(self, dm=None):

    if dm is None: dm = self.dm()

    out = 2*self.ht_obasis.compute_grid_density_dm(
      dm, self.ht_grid.points
    )
    return out

  def Fock_matrix(self, C=None, orthogonalize=False):
    dm = self.dm(C)

    J_kernel = np.tensordot(dm, self.er, axes=([0,1], [0,2]))
    X_kernel = np.tensordot(dm, self.er, axes=([0,1], [0,1]))

    h1 = (self.kin + self.v_ext)
    G = J_kernel * 2. - X_kernel
    F = h1 + G

    if orthogonalize:
      F = self.X.T.dot(F.dot(self.X))

    return F

  def Hamiltonian(self, C=None):
    dm = self.dm(C)

  def get_rho_cube(self, margin=3, resolution=0.1, dm=None):

    if dm is None: dm = self.dm()

    x_min, y_min, z_min = self.molecule.R.min(axis=0) - margin
    x_max, y_max, z_max = self.molecule.R.max(axis=0) + margin

    X, Y, Z = np.meshgrid(
      np.arange(x_min, x_max, resolution),
      np.arange(y_min, y_max, resolution),
      np.arange(z_min, z_max, resolution),
      indexing='ij'
    )

    cube_grid = np.array(zip(X.ravel(), Y.ravel(), Z.ravel()))
    cube_data_list = 2*self.ht_obasis.compute_grid_density_dm(dm, cube_grid)
    cube_data = cube_data_list.reshape(X.shape)

    grid = np.array([
      [self.molecule.N, x_min, y_min, z_min],
      [X.shape[0], resolution, 0, 0],
      [X.shape[1], 0, resolution, 0],
      [X.shape[2], 0, 0, resolution],
    ])

    q = qtk.CUBE()
    q.build(self.molecule, grid, cube_data)

    return q

  def delete_ht_types(self):
    for attr in dir(self):
      if attr.startswith('ht_'):
        delattr(self, attr)
    if type(self.setting['basis_set']) is not str:
      self.setting['basis_set'] = self.setting['basis_set'].filename

  def delete_matrices(self):
    del_list = ['initial_dm', 'initial_mov', 'olp', 'kin', 'na', 'v_ext', 'er']
    for attr in del_list:
      delattr(self, str(attr))

  def matrices(self):
    try:
      return self.initial_dm, self.initial_mov, self.olp, self.kin, self.v_ext, self.er
    except:
      qtk.warning("matrices not found")

class out(GaussianBasisOutput):
  def __init__(self):
    pass
