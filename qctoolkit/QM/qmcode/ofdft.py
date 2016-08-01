import qctoolkit as qtk
import pkgutil
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import periodictable as pt
from libxc import libxc

ps_eggs_loader = pkgutil.find_loader('pyscf')
ps_found = ps_eggs_loader is not None
ht_eggs_loader = pkgutil.find_loader('horton')
ht_found = ht_eggs_loader is not None

if ps_found:
  from pyscf import gto
  import pyscf as ps
else:
  pass
if ht_found:
  from horton import BeckeMolGrid
else:
  pass

dot = np.dot
diag = np.diag
outer = np.outer
sqrt = np.sqrt
eig = np.linalg.eigh

class inp(GaussianBasisInput):
  def __init__(self, molecule, **kwargs):
    if not ht_found:
      qtk.exit("horton module not found.")
    if not ps_found:
      qtk.exit("pyscf module not found.")
    if 'wf_convergence' not in kwargs:
      kwargs['wf_convergence'] = 1e-06
    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()

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

    dv = sqrt(diag(ig))
    dv = dv / sqrt(sum(diag(outer(dv, dv).dot(ovl))))
    dv = dv * sqrt(mol.nelectron)
    dm = outer(dv, dv)

    coord = np.array(np.atleast_2d(molecule.R*1.8897261245650618))
    grid = BeckeMolGrid(coord, molecule.Z.astype(int), molecule.Z)

    self.molecule = molecule
    self.mol = mol
    self.ovl = ovl
    self.kin = kin
    self.ext = ext
    self.vee = vee
    self.nao = nao
    self.initial_guess = copy.deepcopy(dv)
    self.ig = ig
    self.dv = dv
    self.dm = outer(dv, dv)
    self.grid = grid
    self.P = diag(self.dm.dot(ovl))
    self.P_milliken = diag(ig.dot(ovl))

    self.population = [0 for i in range(molecule.N)]
    self.population_lowdin = [0 for i in range(molecule.N)]
    self.population_milliken = [0 for i in range(molecule.N)]
    itr = 0
    for i in range(len(mol._bas)):
      atom = mol.bas_atom(i)
      s = i + itr
      for j in range(2*mol.bas_angular(i)+1):
        self.population[atom] += self.P[s]
        self.population_milliken[atom] += self.P_milliken[s]
        s += 1 
      itr += j

  def getRho(self, coords):
    psi = self.mol.eval_gto("GTOval_sph", coords).T
    rho = np.zeros(len(coords))
    for i in range(len(self.dv)):
      n_i = self.dv[i]
      psi_i = psi[i]
      rho += n_i * psi_i
    rho = rho**2
    return rho

  def getRhoSigma(self, coords):
    psi = self.mol.eval_gto("GTOval_sph", coords).T
    dpsi = self.mol.eval_gto("GTOval_ip_sph", coords, comp=3).T
    rho = np.zeros(len(coords))
    for i in range(len(self.dv)):
      n_i = self.dv[i]
      psi_i = psi[i]
      rho += n_i * psi_i
    rho = rho**2

    drho = np.zeros([len(coords), 3])
    for i in range(len(self.dv)):
      n_i = self.dv[i]
      psi_i = psi[i]
      for j in range(len(self.dv)):
        n_j = self.dv[j]
        dpsi_j = dpsi[j]
        drho += 2 * n_i * n_j * dpsi_j * psi_i[:, np.newaxis]
    sigma = np.sum(drho**2, axis=1)

    return rho, sigma

  def xc(self, xcFlag=1):
    coords = self.grid.points
    rho, sigma = self.getRhoSigma(coords)
    libxc(rho, sigma, len(coords), 1)

  def getCubeGrid(self, margin=7, step=0.2):
    coord_min = np.min(self.mol.atom_coords(), axis=0) - margin
    coord_max = np.max(self.mol.atom_coords(), axis=0) + margin
    steps = np.ceil((coord_max - coord_min)/float(step))
    out = np.zeros([4, 4])
    out[0, 0] = self.molecule.N
    out[0, 1:] = coord_min
    out[1:, 0] = steps
    out[1:, 1:] = diag((coord_max - coord_min) / steps)
    return out

  def dm2cube(self, dv=None, cube_header=None):
    if not dv:
      dv = self.dv
    if not cube_header:
      cube_header = self.getCubeGrid()

    def getSpace(i):
      coord_min = cube_header[0, 1+i]
      step = cube_header[1+i, 0]
      coord_max = coord_min + step*cube_header[1+i, 1+i]
      return np.linspace(coord_min, coord_max, step)

    coord_axes = [getSpace(i) for i in range(3)]
    X, Y, Z = np.meshgrid(*coord_axes, indexing='ij')
    step = X.shape
    X = X.reshape(X.size)
    Y = Y.reshape(Y.size)
    Z = Z.reshape(Z.size)
    coords = np.array([X,Y,Z]).T
    rho = self.getRho(coords)
    rho = rho.reshape(*step)
    cube = qtk.CUBE()
    cube.build(self.molecule, cube_header, rho)

    return cube
  
