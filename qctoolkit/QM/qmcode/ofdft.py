import qctoolkit as qtk
import pkgutil
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import periodictable as pt
from libxc_dict import xc_dict

ps_eggs_loader = pkgutil.find_loader('pyscf')
ps_found = ps_eggs_loader is not None
ht_eggs_loader = pkgutil.find_loader('horton')
ht_found = ht_eggs_loader is not None
xc_eggs_loader = pkgutil.find_loader('libxc')
xc_found = ps_eggs_loader is not None

if ps_found:
  from pyscf import gto
  import pyscf as ps
else:
  qtk.warning("pyscf not found.")
if ht_found:
  from horton import BeckeMolGrid
else:
  qtk.warning("horton not found.")
if xc_found:
  from libxc_exc import libxc_exc
  from libxc_vxc import libxc_vxc
else:
  qtk.warning("libxc not found.")

dot = np.dot
diag = np.diag
outer = np.outer
sqrt = np.sqrt
eig = np.linalg.eigh
td = np.tensordot
trace = np.trace

class inp(GaussianBasisInput):
  def __init__(self, molecule, **kwargs):
    if not ht_found:
      qtk.exit("horton module not found.")
    if not ps_found:
      qtk.exit("pyscf module not found.")
    if not xc_found:
      qtk.exit("libxc not found.")
    if 'wf_convergence' not in kwargs:
      kwargs['wf_convergence'] = 1e-06

    if 'kinetic_functional' not in kwargs:
      kwargs['kinetic_functional'] = 'LLP'

    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()

    if 'dft_setting' not in kwargs:
      kf = self.setting['kinetic_functional']
      if kf == 'LLP':
        self.setting['dft_setting'] = {
          'K': {1.0: 'XC_GGA_K_LLP'},
        }
        
      ss = self.setting['theory']
      sdft = self.setting['dft_setting']
      if ss == 'pbe':
        sdft['X'] = 'XC_GGA_X_PBE'
        sdft['C'] = 'XC_GGA_C_PBE'
      elif ss == 'blyp':
        sdft['X'] = 'XC_GGA_X_B88'
        sdft['C'] = 'XC_GGA_C_LYP'

    dft = self.setting['dft_setting']
    for k, v in dft.iteritems():
      if type(dft[k]) is str:
        dft[k] = xc_dict[v]
      elif type(dft[k]) is dict:
        for key, value in dft[k].iteritems():
          if type(value) is str:
            dft[k][key] = xc_dict[value]

    mol_str = []
    for i in range(molecule.N):
      atm_str = [molecule.type_list[i]]
      for j in range(3):
         atm_str.append(str(molecule.R[i,j]))
      mol_str.append(' '.join(atm_str))
    mol_str = '; '.join(mol_str)

    mol = gto.Mole()
    mol.build(atom=mol_str, basis=self.setting['basis_set'])
    ig = ps.scf.hf.get_init_guess(mol, key = 'atom')
    ovl = mol.intor_symmetric("cint1e_ovlp_sph")
    kin = mol.intor_symmetric("cint1e_kin_sph")
    ext = mol.intor_symmetric("cint1e_nuc_sph")
    # 4D array, can be contracted by numpy.tensordot (td)
    # nao x nao integrated density Coulomb kernel can be computed via
    # int_vee_rho = td(dm, vee, axes=([0,1], [0,1]))
    # the final electron-electron cam be computed via
    # Vee = 0.5 * trace(dm, int_vee_rho)
    nao = mol.nao_nr()
    vee = mol.intor(
            "cint2e_sph", comp=1, hermi=1,
          ).reshape(nao, nao, nao, nao)


    coord = np.array(np.atleast_2d(molecule.R*1.8897261245650618))
    grid = BeckeMolGrid(coord, molecule.Z.astype(int), molecule.Z)

    self.molecule = molecule
    self.mol = mol
    self.ovl = ovl
    self.kin = kin
    self.ext = ext
    self.vee = vee
    self.nao = nao
    self.dv = self.normalize(sqrt(diag(ig)))
    self.initial_guess = copy.deepcopy(self.dv)
    self.ig = ig
    self.grid = grid
    self.update(self.dv)
    self.old_deltadv = None
    self.new_deltadv = None

  def update(self, dv):
    self.dv = dv
    self.dm = outer(dv, dv)
    self.P = diag(self.dm.dot(self.ovl))
    self._phi = None
    self._dphi = None
    self._psi = None
    self._dpsi = None
    self._rho = None
    self._drho = None
    self._sigma = None
    self.population = [0 for i in range(self.molecule.N)]
    mol = self.mol
    itr = 0
    for i in range(len(mol._bas)):
      atom = mol.bas_atom(i)
      s = i + itr
      for j in range(2*mol.bas_angular(i)+1):
        self.population[atom] += self.P[s]
        s += 1 
      itr += j

  def normalize(self, dv):
    dv = dv / sqrt(sum(diag(outer(dv, dv).dot(self.ovl))))
    dv = dv * sqrt(self.mol.nelectron)
    return dv

  def getPhi(self, coords, **kwargs):
    if 'new' in kwargs and kwargs['new']:
      out = self.mol.eval_gto("GTOval_sph", coords).T
    else:
      if self._phi is None:
        self._phi = self.mol.eval_gto("GTOval_sph", coords).T
      out = self._phi
    return out

  def getDphi(self, coords):
    if self._dphi is None:
      self._dphi = self.mol.eval_gto("GTOval_ip_sph", coords, comp=3).T
    return self._dphi

  def getPsi(self, coords=None, **kwargs):
    if coords is None:
      coords = self.grid.points
    if 'new' in kwargs and kwargs['new']:
      phi = self.getPhi(coords, new=True)
      out = np.zeros(len(coords))
      for i in range(len(self.dv)):
        c_i = self.dv[i]
        phi_i = phi[i]
        out += c_i * phi_i
    else:
      if self._psi is None:
        phi = self.getPhi(coords)
        psi = np.zeros(len(coords))
        for i in range(len(self.dv)):
          c_i = self.dv[i]
          phi_i = phi[i]
          psi += c_i * phi_i
        self._psi = psi
      out = self._psi
    return out

  def getDpsi(self, coords=None):
    if self._dpsi is None:
      if coords is None:
        coords = self.grid.points
      dphi = self.getDphi(coords)
      dpsi = np.zeros([len(coords), 3])
      for i in range(len(self.dv)):
        c_i = self.dv[i]
        dphi_i = dphi[i]
        dpsi += c_i * dphi_i
      self._dpsi = dpsi
    return self._dpsi

  def getRho(self, coords=None, **kwargs):
    if coords is None:
      coords = self.grid.points
    if 'new' in kwargs and kwargs['new']:
      out = self.getPsi(coords, new=True)**2
    else:
      if self._rho is None:
        self._rho = self.getPsi(coords)**2
      out = self._rho
    return out

  def getDrho(self, coords=None):
    if self._drho is None:
      if coords is None:
        coords = self.grid.points
      psi = self.getPsi(coords)
      dpsi = self.getDpsi(coords)
      self._drho = 2 * dpsi * psi[:, np.newaxis]
    return self._drho

  def getSigma(self, coords=None):
    if self._sigma is None:
      if coords is None:
        coords = self.grid.points
      drho = self.getDrho(coords)
      self._sigma = np.sum(drho**2, axis=1)
    return self._sigma

  def drhoSigma_dc(self, i, coords = None):
    if coords is None:
      coords = self.grid.points

    phi = self.getPhi(coords)
    dphi = self.getDphi(coords)
    psi = self.getPsi(coords)
    dpsi = self.getDpsi(coords)
    drho = self.getDrho(coords)

    def drho_dc(i):
      return 2 * phi[i] * psi

    def nabla_drho_dc(i):
      term1 = 2 * dphi[i] * psi[:, np.newaxis]
      term2 = 2 * dpsi * phi[i][:, np.newaxis]
      return term1 + term2

    def dsigma_dc(i):
      return np.sum(2 * drho * nabla_drho_dc(i), axis=1)

    return drho_dc(i), dsigma_dc(i)

  def dE_ddv(self, coords = None):
    if coords is None:
      coords = self.grid.points

    rho = self.getRho(coords)
    sigma = self.getSigma(coords)
    dft = self.setting['dft_setting']

    dEk_drho = np.zeros(len(coords))
    dEk_dsigma = np.zeros(len(coords))
    Ek = np.zeros(len(coords))
    for fraction, kf in dft['K'].iteritems():
      dEk_drho_f, dEk_dsigma_f = self.vxc(kf, [rho, sigma], False)
      dEk_drho += fraction * dEk_drho_f
      dEk_dsigma += fraction * dEk_dsigma_f
      Ek += self.exc(kf, [rho, sigma], False)

    dEc_drho, dEc_dsigma = self.vxc(dft['C'], [rho, sigma], False)
    dEx_drho, dEx_dsigma = self.vxc(dft['X'], [rho, sigma], False)
    Ec = self.exc(dft['C'], [rho, sigma], False)
    Ex = self.exc(dft['X'], [rho, sigma], False)
    dE_drho = dEk_drho + dEc_drho + dEx_drho
    dE_dsigma = dEk_dsigma + dEc_dsigma + dEx_dsigma
    epsilon = Ek + Ec + Ex

    dE_kxc = np.zeros(len(self.dv))
    for i in range(len(self.dv)):
      drho_dci, dsigma_dci = self.drhoSigma_dc(i)
      rho_dE_dc = dE_drho * drho_dci + dE_dsigma * dsigma_dci
      E_drho_dc = drho_dci * epsilon
      integrand = rho_dE_dc
      dE_kxc[i] = self.grid.integrate(integrand)

    vee_rho = td(self.dm, self.vee, axes=([0,1], [0,1]))
    dE_ee = 2 * vee_rho.dot(self.dv)
    dE_ext = 2 * self.ext.dot(self.dv)
    out = dE_kxc + dE_ee + dE_ext

    return out

  def E(self, coords = None, **kwargs):
    if coords is None:
      coords = self.grid.points
    rho = self.getRho(coords)
    sigma = self.getSigma(coords)
    dft = self.setting['dft_setting']
    if 'term' not in kwargs:
      e_k = np.zeros(len(coords))
      for fraction, kf in dft['K'].iteritems():
        e_k += fraction * self.exc(kf, [rho, sigma], False)
      K = self.grid.integrate(rho*e_k)
      V = trace(self.dm.dot(self.ext))
      int_vee_rho = td(self.dm, self.vee, axes=([0,1], [0,1]))
      U = trace(self.dm.dot(int_vee_rho))/2.
      intc = self.exc(dft['C'], [rho, sigma], False)
      intx = self.exc(dft['X'], [rho, sigma], False)
      XC = self.grid.integrate(rho * (intc + intx))
      E = K + V + U + XC
    elif kwargs['term'] == 'K':
      e_k = np.zeros(len(coords))
      for fraction, kf in dft['K'].iteritems():
        e_k += fraction * self.exc(kf, [rho, sigma], False)
      E = self.grid.integrate(rho*e_k)
    elif kwargs['term'] == 'vw':
      E = trace(self.dm.dot(self.kin))
    elif kwargs['term'] == 'U':
      int_vee_rho = td(self.dm, self.vee, axes=([0,1], [0,1]))
      E = trace(self.dm.dot(int_vee_rho))/2.
    elif kwargs['term'] == 'V':
      E = trace(self.dm.dot(self.ext))
    elif kwargs['term'] == 'X':
      intx = self.exc(dft['X'], [rho, sigma], False)
      E = self.grid.integrate(rho * intx)
    elif kwargs['term'] == 'C':
      intc = self.exc(dft['C'], [rho, sigma], False)
      E = self.grid.integrate(rho * intc)
    elif kwargs['term'] in xc_dict.values() \
    or kwargs['term'] in xc_dict:
      epsilon = self.exc(kwargs['term'], [rho, sigma])
      E = self.grid.integrate(rho*epsilon)
    return E

  def iterate(self, size=0.005):
    dv = self.dv
    oldE = self.E()
    if self.new_deltadv is None:
      self.old_deltadv = self.dE_ddv()
    else:
      self.old_deltadv = self.new_deltadv
    self.dv = self.normalize(dv - self.old_deltadv*size)
    self.update(self.dv)
    self.new_deltadv = self.dE_ddv()
    newE = self.E()
    print newE, np.dot(self.old_deltadv, self.new_deltadv),\
          np.linalg.norm(self.old_deltadv)

  def libxc_report(self, xc_id, flag):
    for k, v in xc_dict.iteritems():
      if v == xc_id:
        key = k
        qtk.report("libxc_%s" % flag, "xc: %s, id: %d\n" % (key, xc_id))
        break

  def exc(self, xcFlag=1, rhoSigma = None, report=True):
    if type(xcFlag) is int:
      if xcFlag not in xc_dict.values():
        qtk.exit("libxc functional id number %d is not valid" % xcFlag)
      else:
        xc_id = xcFlag
    elif type(xcFlag) is str:
      if xcFlag not in xc_dict:
        qtk.exit("libxc functional id %s is not valid" % xcFlag)
      else:
        xc_id = xc_dict[xcFlag]
    if report:
      self.libxc_report(xc_id, 'exc')
    coords = self.grid.points
    if rhoSigma is None:
      rho = self.getRho(coords)
      sigma = self.getSigma(coords)
    else:
      rho, sigma = rhoSigma
    return libxc_exc(rho, sigma, len(coords), xc_id)

  def vxc(self, xcFlag=1, rhoSigma = None, report=True):
    if type(xcFlag) is int:
      if xcFlag not in xc_dict.values():
        qtk.exit("libxc functional id number %d is not valid" % xcFlag)
      else:
        xc_id = xcFlag
    elif type(xcFlag) is str:
      if xcFlag not in xc_dict:
        qtk.exit("libxc functional id %s is not valid" % xcFlag)
      else:
        xc_id = xc_dict[xcFlag]
    if report:
      self.libxc_report(xc_id, 'vxc')
    coords = self.grid.points
    if rhoSigma is None:
      rho = self.getRho(coords)
      sigma = self.getSigma(coords)
    else:
      rho, sigma = rhoSigma
    return libxc_vxc(rho, sigma, len(coords), xc_id)

  def getCubeGrid(self, margin=7, step=0.2):
    coord_min = np.min(self.mol.atom_coords(), axis=0) - margin
    coord_max = np.max(self.mol.atom_coords(), axis=0) + margin
    steps = np.ceil((coord_max - coord_min)/float(step)) + 1
    out = np.zeros([4, 4])
    out[0, 0] = self.molecule.N
    out[0, 1:] = coord_min
    out[1:, 0] = steps
    out[1:, 1:] = diag((coord_max - coord_min) / (steps - 1))
    return out

  def getCube(self, dv=None, cube_header=None, **kwargs):
    if not dv:
      dv = self.dv
    if not cube_header:
      if 'resolution' not in kwargs:
        res = 0.2
      else:
        res = float(kwargs['resolution'])
      if 'margin' not in kwargs:
        margin = 7.0
      else:
        margin = float(kwargs['resolution'])
      cube_header = self.getCubeGrid(margin, res)

    def getSpace(i):
      coord_min = cube_header[0, 1+i]
      step = cube_header[1+i, 0]
      size = cube_header[1+i, 1+i]
      coord_max = coord_min + (step-1)*size
      return np.linspace(coord_min, coord_max, step)

    coord_axes = [getSpace(i) for i in range(3)]
    X, Y, Z = np.meshgrid(*coord_axes, indexing='ij')
    step = X.shape
    X = X.reshape(X.size)
    Y = Y.reshape(Y.size)
    Z = Z.reshape(Z.size)
    coords = np.array([X,Y,Z]).T
    rho = self.getRho(coords, new=True)
    rho = rho.reshape(*step)
    cube = qtk.CUBE()
    cube.build(self.molecule, cube_header, rho)

    return cube
