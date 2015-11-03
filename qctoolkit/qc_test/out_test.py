#!/usr/bin/python

import qctoolkit as qtk

cpmdOut = qtk.QMOut('data/qmout/cpmd.out')
vaspOut = qtk.QMOut('data/qmout/vasprun.xml', program='vasp')
print cpmdOut.inUnit('eV')
print vaspOut.inUnit('eV')
print cpmdOut
print vaspOut

print cpmdOut - vaspOut

print cpmdOut -'test'
