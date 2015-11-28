#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import gcMatrix
import numpy as np


#nwOut = qtk.QMOut('data/qmout/nwchem/water/anDIR-01_00A.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H2_1g/h.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H_3g/h.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H_2g/h.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H_2a_gaussian_basis/h.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H_1g/h.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H_1a/h.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H2_3g/h2.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H2He_3g/h2he.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H2He_3g-1p/h2he.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H2He_1p/hhe.out', program='nwchem')
nwOut = qtk.QMOut('data/qmout/nwchem/HHe_1p/hhe.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/HHe_2p/hhe.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/HLi_3g/hli.out', program='nwchem')
#nwOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz/h2o.out', program='nwchem')
for b in nwOut.basis:
  print b



#gc.gcint(nwOut.basis, nwOut.R_bohr, nwOut.Z)
#print basisData(nwOut.basis)

#mo = nwOut.mo_vectors
#S = np.array([[1.000000, 0.823248],[0.823248,1.000000]])
#V, U = np.linalg.eig(S)
#X = np.dot(U, np.diag(np.sqrt(1/V)))
#print S
#print U
#print X
#
#
#mox = np.dot(mo, X.T)
#moxS = np.dot(mox, np.dot(S,mox.T))
#moS = np.dot(mo, np.dot(S, mo.T))
#xmoS = np.dot(X, np.dot(moS, X.T))
#
#print np.dot(mox,mox.T)

##print nwOut.Z
aoV = gcMatrix(nwOut.basis, nwOut.R_bohr, nwOut.Z)
moV = np.dot(nwOut.mo_vectors,np.dot(aoV, nwOut.mo_vectors.T))
print aoV
#aoVx = np.dot(X.T, np.dot(aoV, X))
#moVx = np.dot(nwOut.mo_vectors,np.dot(aoVx, nwOut.mo_vectors.T))
#
#Sx = np.dot(X.T, np.dot(S, X))
#moS = np.dot(nwOut.mo_vectors,np.dot(S, nwOut.mo_vectors.T))
#moSx = np.dot(nwOut.mo_vectors,np.dot(Sx, nwOut.mo_vectors.T))
#
#print moVx
##print moV
noc = len([a for a in nwOut.occupation if a > 0])
##print noc
#
print np.dot(np.diag(moV), nwOut.occupation)
##print nwOut.occupation
#print nwOut.basis[10]
#print "aaa"

#print a + 1

#print moS
#print moSx
#print S
#print np.dot(X.T,X)
#print 1/U
#print U
