import qctoolkit as qtk
import glob, re, os, copy, shutil
import numpy as np

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
