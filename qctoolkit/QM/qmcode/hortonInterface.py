import qctoolkit as qtk
import itertools
import pkgutil
eggs_loader = pkgutil.find_loader('horton')
found = eggs_loader is not None
if found:
  try:
    from horton import *
    import horton as ht
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
import qctoolkit.QM.ofdft.libxc_interface as xc_io
import qctoolkit.QM.ofdft.libxc_dict as xc_info

inpPath = os.path.realpath(__file__)
inpPath = os.path.split(inpPath)[0]
xcPath = os.path.join(inpPath, '../ofdft/libxc_fxc.so')
xc_found = os.path.exists(xcPath)
if xc_found:
  from qctoolkit.QM.ofdft.libxc_exc import libxc_exc
  from qctoolkit.QM.ofdft.libxc_vxc import libxc_vxc
  from qctoolkit.QM.ofdft.libxc_fxc import libxc_fxc
  pass
else:
  pass

def is_gga(inp):
  if inp.setting['theory'] in ['pbe', 'blyp']:
    return True
  else:
    return False

def is_hgga(inp):
  if inp.setting['theory'] in ['pbe0', 'b3lyp', 'hse06']:
    return True
  else:
    return False

def is_mgga(inp):
  if inp.setting['theory'] in ['tpss']:
    return True
  else:
    return False
  
def is_mhgga(inp):
  if inp.setting['theory'] in ['m05']:
    return True
  else:
    return False
  

class inp(GaussianBasisInput):
  """
  horton input class. 
  """
  __doc__ = GaussianBasisInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):

    inp.ht = ht

    if 'theory' not in kwargs:
      kwargs['theory'] = 'hf'

    if 'save_c_type' not in kwargs:
      kwargs['save_c_type'] = True

    if not found:
      qtk.exit("horton module not found.")
    if 'wf_convergence' not in kwargs:
      kwargs['wf_convergence'] = 1e-06

    if 'cholesky' not in kwargs:
      kwargs['cholesky'] = True

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
    if self.setting['cholesky']:
      er = obasis.compute_electron_repulsion_cholesky()
    else:
      er = obasis.compute_electron_repulsion()

    exp_alpha = orb_alpha = Orbitals(obasis.nbasis)
    
    # Initial guess
    guess_core_hamiltonian(olp, kin + na, exp_alpha)

    external = {'nn': compute_nucnuc(mol.coordinates, 
                                     mol.pseudo_numbers)}
    theory = self.setting['theory']
    if theory in xc_info.xc_map:
      self.xc = xc_info.xc_map[theory]
    else:
      self.xc = None

    terms = [
       RTwoIndexTerm(kin, 'kin'),
       RTwoIndexTerm(na, 'ne'),
       RDirectTerm(er, 'hartree'),
    ]
    if self.setting['theory'] == 'hf':
      terms.append(RExchangeTerm(er, 'x_hf'))
    elif self.setting['theory'] == 'pbe':
      libxc_terms = [
        RLibXCGGA('x_pbe'),
        RLibXCGGA('c_pbe'),
      ]
      terms.append(RGridGroup(obasis, grid, libxc_terms))
    elif self.setting['theory'] == 'blyp':
      libxc_terms = [
        RLibXCGGA('x_b88'),
        RLibXCGGA('c_lyp'),
      ]
      terms.append(RGridGroup(obasis, grid, libxc_terms))
    elif self.setting['theory'] == 'pbe0':
      hyb_term = RLibXCHybridGGA('xc_pbeh')
      terms.append(RGridGroup(obasis, grid, [hyb_term]))
      terms.append(RExchangeTerm(er, 'x_hf', hyb_term.get_exx_fraction()))
    elif self.setting['theory'] == 'b3lyp':
      hyb_term = RLibXCHybridGGA('xc_b3lyp')
      terms.append(RGridGroup(obasis, grid, [hyb_term]))
      terms.append(RExchangeTerm(er, 'x_hf', hyb_term.get_exx_fraction()))
    elif self.setting['theory'] == 'hse06':
      hyb_term = RLibXCHybridGGA('xc_hse06')
      terms.append(RGridGroup(obasis, grid, [hyb_term]))
      terms.append(RExchangeTerm(er, 'x_hf', hyb_term.get_exx_fraction()))
    elif self.setting['theory'] == 'tpss':
      libxc_terms = [
        RLibXCMGGA('x_tpss'),
        RLibXCMGGA('c_tpss'),
      ]
      terms.append(RGridGroup(obasis, grid, libxc_terms))
    elif self.setting['theory'] == 'm05':
      hyb_term = RLibXCHybridMGGA('xc_m05')
      terms.append(RGridGroup(obasis, grid, [hyb_term]))
      terms.append(RExchangeTerm(er, 'x_hf', hyb_term.get_exx_fraction()))
    ham = REffHam(terms, external)

    occ_model = AufbauOccModel(
      int((sum(self.molecule.Z) - self.molecule.charge)/ 2.)
    )
    occ_model.assign(orb_alpha)

    occ = np.zeros(olp.shape[0])
    N = int(sum(self.molecule.Z) - self.molecule.charge)
    occ[:N/2] = 1
    C = copy.deepcopy(exp_alpha.coeffs.__array__())
    dm = (C * occ).dot(C.T)

    if self.setting['save_c_type']:
      self.ht_mol = mol
      self.ht_grid = grid
      self.grid = grid
      self.ht_external = external
      self.ht_obasis = obasis
      self.ht_occ_model = occ_model
      self.ht_olp = olp
      self.ht_kin = kin
      self.ht_na = na
      self.ht_er = er
      self.ht_exp_alpha = exp_alpha
      self.ht_dm_alpha = exp_alpha.to_dm()
      self.ht_terms = terms
      self.ht_ham = ham
    else:
      self.grid_points = grid.points
      self.grid_weights = grid.weights

    self.olp = olp.__array__()
    self.S = self.olp
    self.kin = kin.__array__()
    self.T = self.kin
    self.nn = external['nn']
    try:
      d, U = np.linalg.eigh(self.olp)
      self.X = U.dot(np.diag(np.sqrt(1/d)).dot(U.T))
    except Exception as err:
      qtk.warning(err)
    self.na = na.__array__()
    self.er = er.__array__()
    self.U = self.er
    self.mov = exp_alpha.coeffs.__array__()
    self.initial_mov = C
    self.initial_dm = dm
    self.occ = occ
    self.occupation = 2*occ
    self.v_ext = self.na
    #self.dm = dm_alpha.__array__()

  def run(self, name=None, **kwargs):

    self.setting.update(kwargs)

    if self.setting['theory'] in ['hf']:
      optimizer = PlainSCFSolver
      opt_arg = [
        self.ht_ham, self.ht_olp, self.ht_occ_model, self.ht_exp_alpha
      ]
    elif is_gga(self) or is_mgga(self):
      optimizer = CDIISSCFSolver
      opt_arg = [
        self.ht_ham, self.ht_olp, self.ht_occ_model, self.ht_dm_alpha
      ]
    elif is_hgga(self) or is_mhgga(self):
      optimizer = EDIIS2SCFSolver
      opt_arg = [
        self.ht_ham, self.ht_olp, self.ht_occ_model, self.ht_dm_alpha
      ]

    scf_solver = optimizer(
      threshold=self.setting['wf_convergence'], 
      maxiter=self.setting['scf_step']
    )
    scf_solver(*opt_arg)

    if self.setting['theory'] is not 'hf':
      fock_alpha = np.zeros(self.ht_olp.shape)
      self.ht_ham.reset(self.ht_dm_alpha)
      self.ht_ham.compute_energy()
      self.ht_ham.compute_fock(fock_alpha)
      self.ht_exp_alpha.from_fock_and_dm(
        fock_alpha, self.ht_dm_alpha, self.ht_olp
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

  def e_kin(self, C=None):
    dm = self.dm(C=C)
    return 2 * np.trace(dm.dot(self.T))

  def e_xc(self, C=None):
    dm = self.dm(C=C)
    exc = self.exc(dm)
    rho = self.getRho(dm)
    return self.ht_grid.integrate(rho * exc)

  def e_coulomb(self, C=None):
    return self.ht_ham.cache['energy_hartree']

  #def e_coulomb(self, C=None):
  #  dm = self.dm(C=C)
  #  J_kernel = np.tensordot(dm, self.U, axes=([0,1], [0,2]))
  #  return 2 * np.trace(dm.dot(J_kernel))

  #def e_xx(self, C=None):
  #  dm = self.dm(C=C)
  #  J_kernel = np.tensordot(dm, self.U, axes=([0,1], [0,1]))
  #  return -np.trace(dm.dot(J_kernel))
    

  def dm(self, C=None):
    if C is None: C = self.mov
    occ = self.occ
    return (C * occ).dot(C.T)

  def Fock_matrix(self, C=None, orthogonalize=False):

    if C is None:
      C = self.mov

    dm = self.dm(C)

    J_kernel = np.tensordot(dm, self.er, axes=([0,1], [0,2]))
    X_kernel = np.tensordot(dm, self.er, axes=([0,1], [0,1]))

    h1 = (self.kin + self.v_ext)
    G = J_kernel * 2. - X_kernel
    F = h1 + G

    if orthogonalize:
      F = C.T.dot(F).dot(C)

    return F

  def cube_grid(self, margin=3, resolution=0.1):
    if margin is None:
      margin = 3
    if resolution is None:
      resolution = 0.1

    x_min, y_min, z_min = self.molecule.R.min(axis=0) - margin
    x_max, y_max, z_max = self.molecule.R.max(axis=0) + margin

    X, Y, Z = np.meshgrid(
      np.arange(x_min, x_max, resolution),
      np.arange(y_min, y_max, resolution),
      np.arange(z_min, z_max, resolution),
      indexing='ij'
    )

    cube_grid = np.array(zip(X.ravel(), Y.ravel(), Z.ravel()))

    grid = np.array([
      [self.molecule.N, x_min, y_min, z_min],
      [X.shape[0], resolution, 0, 0],
      [X.shape[1], 0, resolution, 0],
      [X.shape[2], 0, 0, resolution],
    ])

    return cube_grid, grid, X.shape

  def getPsi(self, C=None, psi_list=None, margin=None, 
             resolution=None, cube=False, get_shape=False, dr=[0,0,0]):
    if C is None: C = self.mov

    mov_back = self.ht_exp_alpha.coeffs.__array__()
    self.ht_exp_alpha._coeffs = C

    exp = self.ht_exp_alpha
    if cube:
      if margin is None: margin = 3
      if resolution is None: resolution = 0.1
    if (margin is None and resolution is None) and not cube:
      dr = np.array(dr)
      assert dr.shape == (3,)
      pts = self.ht_grid.points + dr
    else:
      pts, cube_grid, shape = self.cube_grid(margin, resolution)
    if psi_list is None:
      psi_list = np.array(range(len(self.occ)))
    psi_list = np.array(psi_list)

    psi = self.ht_obasis.compute_grid_orbitals_exp(exp, pts, psi_list)
    
    self.ht_exp_alpha._coeffs = mov_back

    if not get_shape:
      return psi
    else:
      return psi, cube_grid, shape

  def getPhi(self):
    C = self.X.dot(np.eye(self.olp.shape[0]))
    return self.getPsi(C=C)

  def getDPsi(self, C=None, psi_list=None, dr=[0,0,0]):
    if C is None: C = self.mov
    mov_back = self.ht_exp_alpha.coeffs.__array__()
    self.ht_exp_alpha._coeffs = C

    exp = self.ht_exp_alpha

    dr = np.array(dr)
    assert dr.shape == (3,)
    pts = self.ht_grid.points + dr

    if psi_list is None:
      psi_list = np.array(range(len(self.occ)))
    psi_list = np.array(psi_list)

    psi = self.ht_obasis.compute_grid_orb_gradient_exp(exp, pts, psi_list)

    self.ht_exp_alpha._coeffs = mov_back

    return psi

  def getRho(self, dm=None, dr=[0,0,0]):

    if dm is None: dm = self.dm()
    dr = np.array(dr)
    assert dr.shape == (3,)

    out = 2*self.ht_obasis.compute_grid_density_dm(
      dm, self.ht_grid.points + dr
    )
    return out

  def getDRho(self, dm=None, dr=[0,0,0]):

    if dm is None: dm = self.dm()
    dr = np.array(dr)
    assert dr.shape == (3,)

    out = 2*self.ht_obasis.compute_grid_gradient_dm(
      dm, self.ht_grid.points + dr
    )
    return out

  def getSigma(self, dm=None, dr=[0,0,0]):
    drho = self.getDRho(dm, dr=dr)
    return (drho**2).sum(axis=1)

  def _xc_call(self, xc_func, dm, xcFlags, sumup, dr):

    if xcFlags is None:
      xcFlags = self.xc
    #xc_id = xc_io.get_xcid(xcFlag)

    if xcFlags:
      outs = []
      for xcFlag, scale in xcFlags.iteritems():
  
        xc_id = xc_io.get_xcid(xcFlag)
    
        rho = self.getRho(dm, dr=dr)
        sigma = self.getSigma(dm, dr=dr)
        coords = self.ht_grid.points
  
        out = np.stack(xc_func(rho, sigma, len(coords), xc_id))
        if not sumup:
          scale = 1.
        outs.append(scale * out)
        
      if sumup:
        return np.stack(outs).sum(axis=0)
      else:
        return np.stack(outs)

  def exc(self, dm=None, xcFlags=None, sumup=True, dr=[0,0,0]):
    return self._xc_call(libxc_exc, dm, xcFlags, sumup, dr)

  def vxc_raw(self, dm=None, xcFlags=None, sumup=True, dr=[0,0,0]):
    return self._xc_call(libxc_vxc, dm, xcFlags, sumup, dr)

  def vxc(self, dm=None, xcFlags=None, sumup=True, dr=[0,0,0]):
    out = self._xc_call(libxc_vxc, dm, xcFlags, sumup, dr)
    if is_gga(self) or is_hgga(self):
      vrho, vsigma = out[0], out[1]
      vsigma_grad = self.gradFD(self.vxc_raw)[:, 1, :]
      rho_grad = self.gradFD(self.getRho)
      rho_div = self.divFD(self.getRho)

      # functional derivative of GGA
      return vrho - 2 * (rho_div * vsigma + (rho_grad * vsigma_grad).sum(axis=0))

  def fxc(self, dm=None, xcFlags=None, sumup=True, dr=[0,0,0]):
    return self._xc_call(libxc_fxc, dm, xcFlags, sumup, dr)

  def getPsiCube(self, n, C=None, margin=3, resolution=0.1):
    psi, grid, shape = self.getPsi(C, [n], margin, resolution, get_shape=True)
    psi = psi.reshape(shape)

    q = qtk.CUBE()
    q.build(self.molecule, grid, psi)

    return q

  def AOMatrix_grid(self, val_grid):
    phis = self.getPhi()
    out = np.zeros(self.olp.shape)
    for i in range(len(out)):
      for j in range(len(out)):
        out[i, j] = self.ht_grid.integrate(phis[:, i] * val_grid * phis[:, j])
    return out
  def gradFD(self, func, func_kwargs={}, eps=1e-4):
    dx_p = func(dr=[ eps, 0, 0], **func_kwargs)
    dx_n = func(dr=[-eps, 0, 0], **func_kwargs)
    dy_p = func(dr=[ 0, eps, 0], **func_kwargs)
    dy_n = func(dr=[ 0,-eps, 0], **func_kwargs)
    dz_p = func(dr=[ 0, 0, eps], **func_kwargs)
    dz_n = func(dr=[ 0, 0,-eps], **func_kwargs)

    grad_x = (dx_p - dx_n) / (2 * eps)
    grad_y = (dy_p - dy_n) / (2 * eps)
    grad_z = (dz_p - dz_n) / (2 * eps)

    out = np.stack([grad_x, grad_y, grad_z])

    return out

  def divFD(self, func, func_kwargs={}, eps=1e-4):
    dx_p = func(dr=[ eps, 0, 0], **func_kwargs)
    dx_n = func(dr=[-eps, 0, 0], **func_kwargs)
    dy_p = func(dr=[ 0, eps, 0], **func_kwargs)
    dy_n = func(dr=[ 0,-eps, 0], **func_kwargs)
    dz_p = func(dr=[ 0, 0, eps], **func_kwargs)
    dz_n = func(dr=[ 0, 0,-eps], **func_kwargs)
    fr0 =  func(**func_kwargs)

    grad2_x = (dx_p - 2 * fr0 + dx_n) / eps ** 2
    grad2_y = (dy_p - 2 * fr0 + dy_n) / eps ** 2
    grad2_z = (dz_p - 2 * fr0 + dz_n) / eps ** 2

    return grad2_x + grad2_y + grad2_z

  def getRhoCube(self, margin=3, resolution=0.1, dm=None):

    if dm is None: dm = self.dm()

    cube_grid, grid, shape = self.cube_grid(margin, resolution)

    cube_data_list = 2*self.ht_obasis.compute_grid_density_dm(dm, cube_grid)
    cube_data = cube_data_list.reshape(shape)

    q = qtk.CUBE()
    q.build(self.molecule, grid, cube_data)

    return q

  def getSigmaCube(self, margin=3, resolution=0.1, dm=None):

    if dm is None: dm = self.dm()


    cube_grid, grid, shape = self.cube_grid(margin, resolution)

    cube_data_list = 2*self.ht_obasis.compute_grid_gradient_dm(
      dm, cube_grid
    )
    cube_data_list = (cube_data_list ** 2).sum(axis=1)
    cube_data = cube_data_list.reshape(shape)

    q = qtk.CUBE()
    q.build(self.molecule, grid, cube_data)

    return q

  def Minv_matrix(self, Coulomb=True, xc=True):

    eps = self.mo_eigenvalues
    N_occ = int(self.occ.sum())
    i_list = range(N_occ)
    a_list = range(N_occ, self.olp.shape[0])
    ia_list = list(itertools.product(i_list, a_list))

    if Coulomb or xc:
      C = self.mov
      U = self.U
      Umo = np.tensordot(C, U, axes=[[0],[0]])
      Umo = np.tensordot(C, Umo, axes=[[0],[1]])
      Umo = np.tensordot(C, Umo, axes=[[0],[2]])
      Umo = np.tensordot(C, Umo, axes=[[0],[3]])
    
    M = np.zeros([len(ia_list), len(ia_list)])
    for s in range(len(ia_list)):
        i, a = ia_list[s]
        de = eps[a] - eps[i]
        M[s, s] += de

        if Coulomb or xc:
          for t in range(len(ia_list)):
              j, b = ia_list[t]

              if Coulomb:
                M[s, t] += Umo[i,a,j,b]

              if xc:
                pass

    return np.linalg.inv(M), ia_list
            
  def getdRho(self, v_ao, margin=None, resolution=None, 
              cube=False, get_shape=False, Coulomb=True, xc=True):
    """get density response with perturbation v_ao expressed in AO basis"""

    C = self.mov

    v_mo = np.tensordot(C, v_ao, axes=[[0],[0]])
    v_mo = np.tensordot(C, v_mo, axes=[[0],[1]])

    M_inv, ia_list = self.Minv_matrix(Coulomb, xc)

    if cube:
      psi, grid, shape = self.getPsi(C=self.mov, margin=margin, 
          resolution=resolution, cube=cube, get_shape=True)
    else:
      psi = self.getPsi(C=self.mov, margin=margin, 
          resolution=resolution, cube=cube)
    drho = np.zeros(psi[:,0].shape)
    for s in range(len(ia_list)):
        i, a = ia_list[s]
        for t in range(len(ia_list)):
            j, b = ia_list[t]
            drho += M_inv[s, t] * v_mo[j, b] * psi[:, i] * psi[:, a]

    if cube:
      drho = drho.reshape(shape)
      q = qtk.CUBE()
      q.build(self.molecule, grid, drho)
      return q
    else:
      return drho

  def d1E(self, v_ao):
    return np.trace(self.dm().dot(v_ao))

  def d2E(self, v_ao, Coulomb=True, xc=True):
    C = self.mov

    v_mo = np.tensordot(C, v_ao, axes=[[0],[0]])
    v_mo = np.tensordot(C, v_mo, axes=[[0],[1]])

    M_inv, ia_list = self.Minv_matrix(Coulomb, xc)

    dE = 0
    for s in range(len(ia_list)):
        i, a = ia_list[s]
        for t in range(len(ia_list)):
            j, b = ia_list[t]
            dE += M_inv[s, t] * v_mo[j, b] * v_mo[i, a]

    return d2E

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
