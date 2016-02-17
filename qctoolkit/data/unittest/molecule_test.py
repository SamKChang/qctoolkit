import qctoolkit as qtk
import glob, re, os
def test_IO():
  path = os.path.realpath(__file__)
  path = re.sub('[a-zA-Z0-9\._]*$', '', path)
  mol_list = glob.glob(path + 'test_data/molecules/*')
  assert len(mol_list) > 0
  mols = []
  for mol_file in mol_list:
    mols.append(qtk.Molecule(mol_file))
  assert len(mol_list) == len(mols)
  for mol in mols:
    print mol
    print mol.R
    print mol.N
    assert mol.N > 0
