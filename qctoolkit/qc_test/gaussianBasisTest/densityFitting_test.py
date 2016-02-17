#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import basisData
from qctoolkit.QM.gaussianbasis_io import veMatrix
from qctoolkit.QM.gaussianbasis_io import eeMatrix
from qctoolkit.QM.gaussianbasis_io import neMatrix
from qctoolkit.QM.gaussianbasis_io import nnMatrix
from qctoolkit.QM.gaussianbasis_io import vnMatrix
from qctoolkit.QM.gaussianbasis_io import densityFitting as df
from qctoolkit.QM.gaussianbasis_io import densityMatrix as dm
from data import *
import numpy as np
from numpy import tensordot as td


print len(qmOut.basis)
for b in qmOut.basis:
  print b

d = df(qmOut)
NE = neMatrix(qmOut.basis)
C_matrix = 0.5*td(d, NE, axes=(0, 0))
D_matrix = dm(qmOut)
vn = vnMatrix(qmOut.basis, qmOut.R, qmOut.Z)
Evn = np.dot(d, vn)
print "Ext: ",
print Evn
Enn = 0.5*np.trace(np.dot(D_matrix, C_matrix))
print "Coulomb (density-fitting): ",
print Enn
Eee = 0
ee = eeMatrix(qmOut.basis)
occ = [i for i in range(qmOut.n_ao) if qmOut.occupation[i]==2][-1] + 1
out = td(mo, ee, axes=(1,0))
out = td(mo, out, axes=(1,1))
out = td(mo, out, axes=(1,2))
out = td(mo, out, axes=(1,3))
Eee = [out[a,a,b,b] for a in range(occ) for b in range(occ)]
print "Coulomb (exact): ",
print sum(Eee)

print "trace of density matrix"
print np.trace(D_matrix)

print d
print path[0]
