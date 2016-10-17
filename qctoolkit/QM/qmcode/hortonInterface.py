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

    if not found:
      qtk.exit("horton module not found.")
    if 'wf_convergence' not in kwargs:
      kwargs['wf_convergence'] = 1e-06
    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()

    mol = IOData(coordinates=molecule.R, numbers=molecule.Z)
    obasis = get_gobasis(mol.coordinates, mol.numbers,
                         self.setting['basis_set'])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, 
                        mol.pseudo_numbers)

    # Create a linalg factory
    lf = DenseLinalgFactory(obasis.nbasis)
    
    # Compute Gaussian integrals
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, 
                                           mol.pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)

    # Create alpha orbitals
    exp_alpha = lf.create_expansion()
    
    # Initial guess
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    external = {'nn': compute_nucnuc(mol.coordinates, 
                                     mol.pseudo_numbers)}
    libxc_term = RLibXCHybridGGA('xc_b3lyp')
    terms = [
        #RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RGridGroup(obasis, grid, [libxc_term]),
        RExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)

    self.ht_mol = mol
    self.ht_grid = grid
    self.ht_external = external
    self.ht_obasis = obasis
    self.ht_lf = lf
    self.ht_olp = olp
    self.ht_kin = kin
    self.ht_na = na
    self.ht_er = er
    self.ht_exp_alpha = exp_alpha
    self.ht_ham = ham

  def run(self, name=None, **kwargs):
    pass
    
  def write(self, name=None, **kwargs):
    pass

class out(GaussianBasisOutput):
  def __init__(self):
    pass
