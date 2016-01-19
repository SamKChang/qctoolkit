#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import gcMatrix
import numpy as np


#qmOut = qtk.QMOut('data/qmout/nwchem/water/anDIR-01_00A.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_3g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2a_gaussian_basis/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_1g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_1p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2g2p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_aug-cc-pvdz/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_cc-pvdz_no-p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_1a/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_3g/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_3g/h2he.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_3g-1p/h2he.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_1p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_1p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_2p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HLi_3g/hli.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz/h2o.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz_b3lyp/h2o.out', program='nwchem')

#qmOut = qtk.QMOut('data/qmout/gaussian/H_1g/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_1p/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_2g2p/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_2a/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2He_3g-1p/H2He.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2_3g/H2.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2_1s1p/H2.out', program='gaussian')
qmOut = qtk.QMOut('data/qmout/gaussian/H2O_aug-cc-pvdz/H2O.out', program='gaussian')

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
aoV = gcMatrix(qmOut.basis, qmOut.R, qmOut.Z)
moV = np.dot(qmOut.mo_vectors,np.dot(aoV, qmOut.mo_vectors.T))
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
#noc = len([a for a in nwOut.occupation if a > 0])
##print noc
#
print np.dot(np.diag(moV), qmOut.occupation)
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
