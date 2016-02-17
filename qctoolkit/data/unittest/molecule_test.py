import qctoolkit as qtk
import glob, re, os

def setup(**kwargs):
  path = os.path.realpath(__file__)
  path = re.sub('[a-zA-Z0-9\._\-]*$', '', path)
  mols = []
  if 'mol' not in kwargs:
    mol_list = glob.glob(path + 'test_data/molecules/*')
  else:
    mol_list = glob.glob(path + 'test_data/molecules/' + kwargs['mol'])
  print path
  print mol_list
  assert len(mol_list) > 0
  for mol_file in mol_list:
    mols.append(qtk.Molecule(mol_file))
  assert len(mol_list) == len(mols)

  return mols

def test_IO():
  # test for read and __add__ function
  mols = setup()
  for mol in mols:
    assert mol.N > 0
    assert mol.R.shape == (mol.N, 3)
  new = mols[0] + mols[1]
  assert new.N == mols[0].N + mols[1].N

def test_h2o_oh_properties():
  mol = setup(mol='h2o-oh.xyz')[0]
  assert mol.N == 5
  assert mol.charge == -1
  assert mol.have_bond('H', 'O')
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

  mol.find_bonds(charge_saturation=1)
  molA = mol.segments[0]
  molB = mol.segments[1]
  assert molA.charge == 0
  assert molB.charge == 1
  assert molA.getValenceElectrons() == 8
  assert molB.getValenceElectrons() == 6

