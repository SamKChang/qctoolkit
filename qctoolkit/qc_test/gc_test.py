#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import gcMatrix
import numpy as np


nwOut = qtk.QMOut('data/qmout/water_nwchem/anDIR-01_00A.out', program='nwchem')
#gc.gcint(nwOut.basis, nwOut.R_bohr, nwOut.Z)
#print basisData(nwOut.basis)
aoV = gcMatrix(nwOut.basis, nwOut.R_bohr, nwOut.Z)
moV = np.dot(nwOut.mo_coefficients,np.dot(aoV, nwOut.mo_coefficients.T))
print aoV
#print sum(np.diag(moV)[:5])
#print nwOut.occupation
#print nwOut.basis[10]
#print "aaa"

#print a + 1


