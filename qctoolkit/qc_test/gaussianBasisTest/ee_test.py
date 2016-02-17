#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import basisData
from qctoolkit.QM.gaussianbasis_io import veMatrix
from qctoolkit.QM.gaussianbasis_io import eeMatrix
import numpy as np
from numpy import tensordot as td
from data import *

#print len(qmOut.basis)
#for b in qmOut.basis:
#  print b

ee = eeMatrix(qmOut.basis)
ex = np.swapaxes(ee, 1,2)
fock = 2*ee - ex

# NOTE: order of contraction matters!
# 0,1,2,3 is correct while 3,2,1,0 is wrong!
out = td(mo, fock, axes=(1,0))
out = td(mo, out, axes=(1,1))
out = td(mo, out, axes=(1,2))
out = td(mo, out, axes=(1,3))
Eee = [out[a,a,b,b] for a in range(occ) for b in range(occ)]
print "tensordot: ",
print sum(Eee)
print path[0]
