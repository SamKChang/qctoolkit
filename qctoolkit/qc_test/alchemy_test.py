#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import gcMatrix
import numpy as np


qmOut = qtk.QMOut('data/qmout/nwchem/alchemy/l_0_basis/l_0.1.out', program='nwchem')
print qmOut.n_ao
print qmOut

#aoV = gcMatrix(qmOut.basis, qmOut.R, qmOut.Z)
#moV = np.dot(qmOut.mo_vectors,np.dot(aoV, qmOut.mo_vectors.T))
#print aoV
#print np.dot(np.diag(moV), qmOut.occupation)
