#!/usr/bin/python

import qctoolkit as qtk

cpmdOut = qtk.QMOut('data/qmout/cpmd.out')
vaspOut = qtk.QMOut('data/qmout/vasp/vasprun.xml', program='vasp')
#nwOut = qtk.QMOut('data/qmout/nwchem/water/anDIR-01_00A.out', program='nwchem')
nwOut = qtk.QMOut('data/qmout/nwchem/A_nwchem/A.out', program='nwchem')
gsOut = qtk.QMOut('data/qmout/gaussian/H2O_aug-cc-pvdz/H2O.out', program='gaussian')
print gsOut.n_mo, gsOut.n_ao
print len(gsOut.basis)
#print cpmdOut.inUnit('eV')
#print vaspOut.inUnit('eV')
#print cpmdOut
#print vaspOut
#
#print cpmdOut - vaspOut
#
#print cpmdOut -'test'

#print nwOut
#print nwOut.n_basis
#print nwOut.mo
#print nwOut.R
#print nwOut.R_bohr
#print nwOut.type_list
#print nwOut.basis_atoms
#print nwOut.basis_type
#print nwOut.basis_exponents
#print nwOut.basis_coefficients
#print nwOut.ao_keys
#print nwOut.basis
#print len(nwOut.basis)
