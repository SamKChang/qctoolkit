#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import veMatrix
from qctoolkit.QM.atomicbasis_io import eeMatrix
from qctoolkit.QM.atomicbasis_io import neMatrix
from qctoolkit.QM.atomicbasis_io import nnMatrix
import numpy as np
from numpy import tensordot as td
from data import *

#print len(qmOut.basis)
#for b in qmOut.basis:
#  print b

nn = nnMatrix(qmOut.basis)
print nn
print nn.shape
