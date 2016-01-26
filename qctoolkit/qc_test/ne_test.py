#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import basisData
from qctoolkit.QM.atomicbasis_io import veMatrix
from qctoolkit.QM.atomicbasis_io import eeMatrix
from qctoolkit.QM.atomicbasis_io import neMatrix
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
#qmOut = qtk.QMOut('data/qmout/nwchem/H2_1g/h.out', program='nwchem')
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
qmOut = qtk.QMOut('data/qmout/gaussian/H2O_aug-cc-pvdz/H2O.out', program='gaussian')

print len(qmOut.basis)
for b in qmOut.basis:
  print b

new_basis = qmOut.basis[0:-1]
for b in new_basis:
  print b
print len(new_basis)


ne = neMatrix(qmOut.basis)
print ne
print ne.shape
mo = qmOut.mo_vectors
mo_matrix = td(
  np.atleast_2d(mo[0,:]), np.atleast_2d(mo[0,:]), 
  axes=(1,1))

print qmOut.occupation
