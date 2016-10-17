import qctoolkit as qtk
import pkgutil
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
from libxc_dict import xc_dict
import libxc_interface as xcio
import copy
import numpy as np
import grid_points as gp
import cube as cube
import os
import evaluation as evaluate
import scipy.optimize as opt

ps_eggs_loader = pkgutil.find_loader('pyscf')
ps_found = ps_eggs_loader is not None
ht_eggs_loader = pkgutil.find_loader('horton')
ht_found = ht_eggs_loader is not None
selfPath = os.path.realpath(__file__)
selfPath = os.path.split(selfPath)[0]
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

dot = np.dot
diag = np.diag
outer = np.outer
sqrt = np.sqrt

class inp(GaussianBasisInput):
  def __init__(self, molecule, **kwargs):
    if not ht_found:
      qtk.exit("horton module not found.")
    if not ps_found:
      qtk.exit("pyscf module not found.")
    if not xc_found:
      print xcpath
      print xc_found
      qtk.exit("libxc not found.")
    if 'wf_convergence' not in kwargs:
      kwargs['wf_convergence'] = 1e-06

    if 'kinetic_functional' not in kwargs:
      kwargs['kinetic_functional'] = 'LLP'

    if 'aufbau' in kwargs and kwargs['aufbau']:
      self.aufbau = True
      self.orbitals = []
    else:
      self.aufbau = False

    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()

    if 'dft_setting' not in kwargs:
      if not self.aufbau:
        kf = self.setting['kinetic_functional']
        if kf == 'LLP':
          self.setting['dft_setting'] = {
            'K': {1.0: 'XC_GGA_K_LLP'},
          }
      else:
        self.setting['dft_setting'] = {
          'K': {1.0: 'XC_GGA_K_VW'}
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
    self.old_E = None
    self.dv_list = []
    self.psi_list = []
    self.dpsi_list = []

    self.optimizers = {
      'tnc': opt.fmin_tnc,
      'ncg': opt.fmin_ncg,
      'cg': opt.fmin_cg,
      'bfgs': opt.fmin_bfgs,
      'l_bfgs_b': opt.fmin_l_bfgs_b,
      'simplex': opt.fmin,
    }
    self.optimizer_settings = {
      'tnc': {'xtol': 0.0, 'pgtol': 0.0, 'maxfun': 1000},
      'simplex': {'xtol': 1E-10, 'ftol': 1E-10, 'maxfun': 1000},
    }

    self.orbitals = []

  def initialize(self):
    self.update(self.initial_guess)

  def update(self, dv):
    self.dv = self.normalize(dv)
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

  def normalize(self, dv=None):
    update = False
    if dv is None:
      update = True
      dv = self.dv
    if self.aufbau:
      N = 2
    else:
      N = self.mol.nelectron
    dv = dv / sqrt(sum(diag(outer(dv, dv).dot(self.ovl))))
    dv = dv * sqrt(N)
    if update:
      self.update(dv)
    return dv

  def E(self, coords = None, **kwargs):
    self.normalize()
    return evaluate.E(self, coords, **kwargs)

  def E_dv(self, dv):
    dv = self.normalize(dv)
    return evaluate.E_dv(self, dv)

  def dE_ddv(self, coords = None):
    return evaluate.dE_ddv(self, coords = None)

  def new_E(self, dv):
    return evaluate.eval_E(self, dv)

  def new_dE(self, dv):
    return evaluate.eval_dE_ddv(self, dv)

  def optimize(self, **kwargs):
    f = self.new_E
    df = self.new_dE
    if self.old_E is None:
      self.old_E = self.E()
    if 'method' in kwargs:
      method = kwargs['method']
    else:
      method = 'tnc'
    if 'tolerance' not in kwargs:
      tolerance = 1E-8
    else:
      tolerance = kwargs['tolerance']
    if 'setting' not in kwargs:
      if method in self.optimizer_settings:
        setting = self.optimizer_settings[method]
      else:
        setting = {}
    else:
      setting = kwargs['setting']
    if 'finit_difference' in kwargs and kwargs['finit_difference']:
      pass
    else:
      if method != 'simplex':
        setting['fprime'] = df

    err, itr = 1, 1
    minE = self.optimizers[method]
    while abs(abs(err) - tolerance) > 1E-8:
      print itr, err, self.old_E
      x0 = copy.deepcopy(self.dv)
      dv_final = minE(f, x0, **setting)
      err = self.E() - self.old_E
      self.old_E += err
      itr += 1
    self.update(self.dv)

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
    norm1 = np.linalg.norm(self.old_deltadv)
    norm2 = np.linalg.norm(self.new_deltadv)
    angle = np.dot(self.old_deltadv, self.new_deltadv) / (norm1 * norm2)
    angle = np.arccos(angle) * 180 / np.pi
    msg = "E=% 13.6E, grad_angle=% 8.3f, |grad|=%8.3f\n" % \
          (newE, angle, norm1)
    qtk.progress("ofdft", msg)

  def CUBE(self, dv=None, cube_header=None, **kwargs):
    return cube.getCube(self, dv, cube_header, **kwargs)
