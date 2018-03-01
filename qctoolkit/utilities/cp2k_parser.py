import qctoolkit as qtk
import numpy as np
import re

def read_ao_file(path):
  with open(path) as f:
    #data = f.readlines()
    data = f.read().splitlines()

  keys = filter(lambda x: 'matrix' in x.lower(), data)

  out = {}
  for i, k in enumerate(keys):
    ind_start = data.index(k)
    ind_end = data.index(keys[i+1]) if i+1 != len(keys) else -1

    key = k.replace('\n', '').lower()[1:]

    matrix_str = filter(lambda x: len(x) > 2, data[ind_start+1:ind_end])
    inds_parts = filter(lambda x: '              ' in x, matrix_str)

    matrix_parts = []
    for j, s in enumerate(inds_parts):
      part_start = matrix_str.index(s) + 1
      if s != inds_parts[-1]:
        part_end = matrix_str.index(inds_parts[j+1])
      else:
        part_end = None
      matrix_part = matrix_str[part_start:part_end]
      part_data = np.array(
        [terms.split()[4:] for terms in matrix_part]
      ).astype(float)
      matrix_parts.append(part_data)

      if k == keys[-1] and s == inds_parts[-1]:
        out['map'] = ["-".join(term.split()[2:4]) for term in matrix_part]

    out[key] = np.hstack(matrix_parts)

  return out

def read_mo_file(path):
  pass

def read_basis_file(path):
  pass
