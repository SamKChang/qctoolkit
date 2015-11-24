#!/usr/bin/python

import qctoolkit as qtk

cpmdOut = qtk.QMOut('data/qmout/cpmd.out')
vaspOut = qtk.QMOut('data/qmout/vasprun.xml', program='vasp')
#nwOut = qtk.QMOut('data/qmout/water_nwchem/anDIR-01_00A.out', program='nwchem')
nwOut = qtk.QMOut('A/A.out', program='nwchem')
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
print nwOut.R
print nwOut.R_bohr
print nwOut.type_list
print nwOut.basis_atoms
print nwOut.shell_list
print nwOut.basis_exponents
print nwOut.basis_coefficients
print nwOut.ao_keys
