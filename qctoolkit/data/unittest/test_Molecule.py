#import qctoolkit as qtk
#import glob, re, os, copy
#import numpy as np
#
#def setup(**kwargs):
#  path = os.path.realpath(__file__)
#  path = re.sub('[a-zA-Z0-9\._\-]*$', '', path)
#  mols = []
#  if 'mol' not in kwargs:
#    mol_list = glob.glob(path + 'test_data/molecules/*')
#  else:
#    mol_list = glob.glob(path + 'test_data/molecules/' + kwargs['mol'])
#  print path
#  print mol_list
#  assert len(mol_list) > 0
#  for mol_file in mol_list:
#    mols.append(qtk.Molecule(mol_file))
#  assert len(mol_list) == len(mols)
#
#  return mols
from setup_test import *

def test_IO():
  # test for read and __add__ function
  mols = setup()
  for mol in mols:
    print mol.N
    print mol.R
    print mol.Z
    print mol.type_list
    print mol
    assert mol.N > 0
    assert mol.R.shape == (mol.N, 3)
  try:
    new = mols[0] + mols[1]
    assert new.N == mols[0].N + mols[1].N
    new.write()
    new.write_xyz()
    new.write_pdb()
    del mols
    del new
  except RuntimeError:
    pass

def test_h2o_oh_segments():
  mol = setup(mol='h2o-oh.xyz')[0]
  assert mol.stoichiometry() == 'H3O2'
  assert mol.N == 5
  assert mol.charge == -1
  assert mol.haveBond('H', 'O')
  assert len(mol.bonds) == 3
  assert mol.bond_types['H-O'] == 3
  molA = mol.segments[0]
  molB = mol.segments[1]
  assert molA.N == 3
  assert molB.N == 2
  assert molA.charge == 0
  assert molB.charge == -1
  assert molA.getValenceElectrons() == 8
  assert molB.getValenceElectrons() == 8

  mol.findBonds(charge_saturation=1)
  molA = mol.segments[0]
  molB = mol.segments[1]
  assert molA.charge == 0
  assert molB.charge == 1
  assert molA.getValenceElectrons() == 8
  assert molB.getValenceElectrons() == 6
  molB.setCharge()
  assert molB.charge == -1
  molB.setCharge(charge_saturation=3)
  assert molB.charge == 3
  del mol
  del molA
  del molB

def test_h2o_center():
  mol = setup(mol='h2o.xyz')[0]
  new1 = copy.deepcopy(mol)
  new2 = copy.deepcopy(mol)
  c = new1.getCenter()
  new1.center(c)
  new2.center(c)
  new2.shift([2,0,0])
  new3 = new1 + new2
  for i in range(2):
    assert new3.distance(i, i+3) == 2
  del mol
  del new1
  del new2
  del new3

def test_h2o_centerOfMass():
  mol = setup(mol='h2o.xyz')[0]
  new1 = copy.deepcopy(mol)
  new2 = copy.deepcopy(mol)
  c = new1.getCenterOfMass()
  new1.center(c)
  new2.center(c)
  new2.shift([2,0,0])
  new3 = new1 + new2
  for i in range(2):
    assert new3.distance(i, i+3) == 2
  del mol
  del new1
  del new2
  del new3

def test_h2o_centerOfCharge():
  mol = setup(mol='h2o.xyz')[0]
  new1 = copy.deepcopy(mol)
  new2 = copy.deepcopy(mol)
  c = new1.getCenterOfCharge()
  new1.center(c)
  new2.center(c)
  new2.shift([2,0,0])
  new3 = new1 + new2
  for i in range(2):
    assert new3.distance(i, i+3) == 2
  del mol
  del new1
  del new2
  del new3

def test_h2o_edit():
  mol = setup(mol='h2o.xyz')[0]
  mol.align([2,1,0])
  crd = np.array([0, 0.87349, 0.28382])
  mol.addAtoms('H', crd)
  assert mol.N == 4
  assert np.linalg.norm(mol.R[3]) == np.linalg.norm(crd)
  mol.addAtoms(['H', 'H'], [crd+1, crd+2])
  assert mol.N == 6
  mol.setAtoms(2, element='S')
  mol.setAtoms([0,1], element='Li')
  assert mol.Z[0] == 3
  assert mol.Z[1] == 3
  assert mol.Z[2] == 16
  mol.setAtoms(2, string='test')
  mol.removeAtoms([3,5,4])
  print mol.N
  assert mol.N == 3
  mol.isolateAtoms([1,2])
  assert mol.N == 2
  assert mol.Z[0] == 3
  assert mol.Z[1] == 0 # artifitial atom with Z<=0

def test_h2o_periodic_crystal():
  mol = setup(mol='h2o.xyz')[0]
  v = mol.R[1] - mol.R[2]
  mol.center(mol.R[2])
  mol.align(v)
  v_new = mol.R[1] - mol.R[2]
  v_norm = np.linalg.norm(v_new)

  assert abs(np.linalg.norm(np.array(v)) - v_norm) < 1E-7
  assert abs(mol.R[1,0] - v_norm) < 1E-7
  assert mol.R[1,2] < 1E-7

  for i in range(1,6):
    mol_new = mol.copy()
    mol_new.stretch(1,[2,1],i)

    assert abs(mol_new.R[1,0] - (v_norm + i)) < 1E-7

def test_crystal():
  mol = setup(mol='periodic_algaas.xyz')[0]
  assert mol.celldm
  assert mol.scale
