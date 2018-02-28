import qctoolkit as qtk
import numpy as np

def read_ao_file(path):
  with open(path) as f:
    data = f.readlines()

  keys = filter(lambda x: 'matrix' in x.lower(), data)

  out = {}
  for i, k in enumerate(keys):
    ind_start = data.index(k)
    ind_end = data.index(keys[i+1]) if i+1 != len(keys) else -1

    key = k.replace('\n', '').lower()[1:]

    matrix_str = filter(lambda x: len(x) > 2, data[ind_start:ind_end])
    out[key] = matrix_str

  return out

def read_mo_file(path):
  pass

def read_basis_file(path):
  pass
