#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import basisData
from qctoolkit.QM.gaussianbasis_io import veMatrix
from qctoolkit.QM.gaussianbasis_io import keMatrix
from qctoolkit.QM.gaussianbasis_io import densityMatrix as dm
import numpy as np
from numpy import tensordot as td
from data import *

print len(qmOut.basis)
for b in qmOut.basis:
  print b

ke = keMatrix(qmOut.basis)
print "kinetic energy matrix"
print ke

D = dm(qmOut)
print "density matrix"
print D
print np.trace(np.dot(D, ke))

print path[0]
