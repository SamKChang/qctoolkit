#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import gcMatrix
from qctoolkit.QM.atomicbasis_io import eeMatrix
import numpy as np


#qmOut = qtk.QMOut('data/qmout/nwchem/H_1g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_1a/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_3g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_1p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_1g1p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2g1p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_3g1p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2g2p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_2a_gaussian_basis/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_aug-cc-pvdz/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H_cc-pvdz_no-p/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g/h.out', program='nwchem')
qmOut = qtk.QMOut('data/qmout/nwchem/H2_1a/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_3g/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_3g/h2he.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_3g-1p/h2he.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_1p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_1p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_2p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HLi_3g/hli.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz/h2o.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz_b3lyp/h2o.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/water/anDIR-01_00A.out', program='nwchem')

#qmOut = qtk.QMOut('data/qmout/gaussian/H_1g/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_1p/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_2g2p/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_2a/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2He_3g-1p/H2He.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2_3g/H2.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2_1s1p/H2.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2O_aug-cc-pvdz/H2O.out', program='gaussian')

print len(qmOut.basis)
for b in qmOut.basis:
  print b

ee = eeMatrix(qmOut.basis)
print ee
#for i in range(10):
#  for j in range(10):
#    for k in range(10):
#      for l in range(10):
#        print i, j, k, l, ee[i,j,k,l]
print ee.shape

occ = [i for i in range(qmOut.n_ao) if qmOut.occupation[i]==2][-1] + 1

mo = qmOut.mo_vectors
Eee = 0
for i in range(ee.shape[0]):
  for j in range(ee.shape[1]):
    for k in range(ee.shape[2]):
      for l in range(ee.shape[3]):
#        print "(%d%d|%d%d):" % (i,j,k,l),
#        print ee[i,j,k,l]
        for a in range(occ):
          for b in range(occ):
#            print mo[a,i]
#            print mo[a,j]
#            print mo[b,k]
#            print mo[b,l]
#            print mo[a,i]*mo[a,j]*mo[b,k]*mo[b,l]*ee[i,j,k,l]
            Eee += mo[a,i]*mo[a,j]*mo[b,k]*mo[b,l]*ee[i,j,k,l]
#        Eee += ee[i,j,k,l] * qmOut.mo_vectors
print qmOut.mo_vectors

print qmOut.n_ao

print Eee
print qmOut.occupation


#print "(00|00):",
#print ee[0,0,0,0]
#print "(11|11):",
#print ee[1,1,1,1]
#print "(00|01):",
#print ee[0,0,0,1]
#print "(00|10):",
#print ee[0,0,1,0]
#print "(01|00):",
#print ee[0,1,0,0]
#print "(10|00):",
#print ee[1,0,0,0]
#print "(00|11):",
#print ee[0,0,1,1]
#print "(11|00):",
#print ee[1,1,0,0]
#print "(10|01):",
#print ee[1,0,0,1]
#print "(10|01):",
#print ee[0,1,0,1]
