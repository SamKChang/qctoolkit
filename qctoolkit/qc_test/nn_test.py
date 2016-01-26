#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import veMatrix
from qctoolkit.QM.atomicbasis_io import eeMatrix
from qctoolkit.QM.atomicbasis_io import neMatrix
from qctoolkit.QM.atomicbasis_io import nnMatrix
import numpy as np
from numpy import tensordot as td


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
qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g-4e/h.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1a/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_3g/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_aug-cc-pvdz/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2-y_aug-cc-pvdz/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2-yz_aug-cc-pvdz/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-z/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1a/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-012/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-013/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-023/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-123/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g-asym-01/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g-asym-02/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g-asym-12/h2.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_aug-cc-pvdz-2e/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_aug-cc-pvdz-4e/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-asym-012/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-asym-012-2e/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-asym-013/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-asym-023/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-asym-123/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H3_1g-asym/h3.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H4_1g/h4.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_1p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_1p-y/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HHe_2p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/HLi_3g/hli.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_1p/hhe.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_3g/h2he.out', program='nwchem')
#qmOut = qtk.QMOut('data/qmout/nwchem/H2He_3g-1p/h2he.out', program='nwchem')
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

#print len(qmOut.basis)
#for b in qmOut.basis:
#  print b

nn = nnMatrix(qmOut.basis)
print nn
print nn.shape
#ex = np.swapaxes(ee, 1,2)
#fock = 2*ee - ex
#occ = [i for i in range(qmOut.n_ao) if qmOut.occupation[i]==2][-1] + 1
#mo = qmOut.mo_vectors
#Eee = 0
#for a in range(occ):
#  for b in range(occ):
#    for i in range(ee.shape[0]):
#      for j in range(ee.shape[1]):
#        for k in range(ee.shape[2]):
#          for l in range(ee.shape[3]):
#            Eee += mo[a,i]*mo[a,j]*mo[b,k]*mo[b,l]*fock[i,j,k,l]
#print "direct for loop: ",
#print Eee
#
## NOTE: order of contraction matters!
## 0,1,2,3 is correct while 3,2,1,0 is wrong!
#out = td(mo, fock, axes=(1,0))
#out = td(mo, out, axes=(1,1))
#out = td(mo, out, axes=(1,2))
#out = td(mo, out, axes=(1,3))
#Eee = [out[a,a,b,b] for a in range(occ) for b in range(occ)]
#print "tensordot: ",
#print sum(Eee)
