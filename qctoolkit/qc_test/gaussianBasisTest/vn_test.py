#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import basisData
from qctoolkit.QM.gaussianbasis_io import vnMatrix
import numpy as np
from data import *

#print len(qmOut.basis)
#for b in qmOut.basis:
#  print b

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
aoV = vnMatrix(qmOut.basis, qmOut.R, qmOut.Z)
print aoV
print aoV.shape
#aoVx = np.dot(X.T, np.dot(aoV, X))
#moVx = np.dot(nwOut.mo_vectors,np.dot(aoVx, nwOut.mo_vectors.T))
#
#Sx = np.dot(X.T, np.dot(S, X))
#moS = np.dot(nwOut.mo_vectors,np.dot(S, nwOut.mo_vectors.T))
#moSx = np.dot(nwOut.mo_vectors,np.dot(Sx, nwOut.mo_vectors.T))
#
#print moVx
##print moV
#noc = len([a for a in nwOut.occupation if a > 0])
##print noc
#
#for i in range(qmOut.n_ao):
#  for j in range(i, qmOut.n_ao):
#    print aoV[i,j]
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
print path[0]
