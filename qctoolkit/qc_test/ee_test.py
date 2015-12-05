#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import gcMatrix
from qctoolkit.QM.atomicbasis_io import eeMatrix
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
qmOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz/h2o.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2O_aug-cc-pvdz_b3lyp/h2o.out', program='nwchem')

#qmOut = qtk.QMOut('data/qmout/gaussian/H_1g/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_1p/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_2g2p/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H_2a/H.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2He_3g-1p/H2He.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2_3g/H2.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2_1s1p/H2.out', program='gaussian')
#qmOut = qtk.QMOut('data/qmout/gaussian/H2O_aug-cc-pvdz/H2O.out', program='gaussian')

#print len(qmOut.basis)
#for b in qmOut.basis:
#  print b

ee = eeMatrix(qmOut.basis)
print ee
#for i in range(10):
#  for j in range(10):
#    for k in range(10):
#      for l in range(10):
#        print i, j, k, l, ee[i,j,k,l]
print ee.shape
