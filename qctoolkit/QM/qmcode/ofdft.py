import qctoolkit as qtk
from pyscf import gto
import pyscf as ps
import pkgutil
eggs_loader = pkgutil.find_loader('horton')
found = eggs_loader is not None
if found:
  from horton import BeckeMolGrid
else:
  pass
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import periodictable as pt

from horton import *

dot = np.dot
diag = np.diag
outer = np.outer
sqrt = np.sqrt
eig = np.linalg.eigh

class inp(GaussianBasisInput):
  def __init__(self, molecule, **kwargs):
    if not found:
      qtk.exit("horton module not found.")
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

    # double check
    ht_mol = IOData(coordinates=molecule.R*1.8897261245650618,
                    numbers=molecule.Z)
    obasis = get_gobasis(ht_mol.coordinates, ht_mol.numbers,
                         self.setting['basis_set'])
    lf = DenseLinalgFactory(obasis.nbasis)
    ht_ovl = obasis.compute_overlap(lf)
    ht_kin = obasis.compute_kinetic(lf)
    ht_rep = obasis.compute_electron_repulsion(lf)
    na = obasis.compute_nuclear_attraction(
           ht_mol.coordinates, ht_mol.pseudo_numbers, lf
         )
    two = obasis.compute_electron_repulsion(lf)
    print ovl.shape
    print ht_kin._array.shape
    print np.sum(abs(kin - ht_kin._array))
    print np.sum(abs(ovl - ht_ovl._array))
    print np.sum(abs(ext - na._array))
    print np.sum(abs(ext))
    print two._array[0,0,0,0]
    print vee[0,0,0,0]
    print np.sum(two._array - vee)

    s, U = eig(ovl)
    sqrS = U.dot(diag(sqrt(s))).dot(U.T)
    X = U.dot(diag(1/sqrt(s))).dot(U.T)
    P = sqrS.dot(ig.dot(sqrS.T))
    dv_orth = sqrt(diag(P))
    dv = X.dot(dv_orth) 

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
    self.dv = dv
    self.dm = outer(dv, dv)
    self.grid = grid
    
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

  def dm2cube(self, dm=None, cube_header=None):
    if not dm:
      dm = self.dm
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
    psi = self.mol.eval_gto("GTOval_sph", coords).T

    rho = np.zeros(X.size)
    for i in range(len(dm)):
      psi_i = psi[i]
      for j in range(len(dm)):
        n_ij = dm[i, j]
        psi_j = psi[j]
       
        rho += n_ij * (psi_i * psi_j)
    rho = rho.reshape(*step)
    cube = qtk.CUBE()
    cube.build(self.molecule, cube_header, rho)

    return cube
  
