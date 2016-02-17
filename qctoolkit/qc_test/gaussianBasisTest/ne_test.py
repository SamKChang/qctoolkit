#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import basisData
from qctoolkit.QM.gaussianbasis_io import veMatrix
from qctoolkit.QM.gaussianbasis_io import eeMatrix
from qctoolkit.QM.gaussianbasis_io import neMatrix
import numpy as np
from numpy import tensordot as td
from data import *

print len(qmOut.basis)
for b in qmOut.basis:
  print b

new_basis = qmOut.basis[0:-1]
for b in new_basis:
  print b
print len(new_basis)

ne = neMatrix(qmOut.basis)
print ne
print ne.shape
mo_matrix = td(
  np.atleast_2d(mo[0,:]), np.atleast_2d(mo[0,:]), 
  axes=(1,1))
print path[0]
