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
  nwchem input class. 
  """
  __doc__ = GaussianBasisInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):

    if 'theory' not in kwargs:
      kwargs['theory'] = 'hf'

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

    occ_model = AufbauOccModel(self.ve() / 2)

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

    self.olp = olp.__array__()
    self.kin = kin.__array__()
    self.na = na.__array__()
    self.er = er.__array__()
    self.mov = exp_alpha.coeffs.__array__()
    #self.dm = dm_alpha.__array__()

  def run(self, name=None, **kwargs):
    scf_solver = PlainSCFSolver(self.setting['wf_convergence'])
    #scf_solver = CDIISSCFSolver(1e-6)
    scf_solver(
      self.ht_ham, 
      self.ht_olp, 
      self.ht_occ_model, 
      self.ht_exp_alpha
    )
    
  def write(self, name=None, **kwargs):
    pass

class out(GaussianBasisOutput):
  def __init__(self):
    pass
