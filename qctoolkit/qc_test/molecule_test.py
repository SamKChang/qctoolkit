#!/usr/bin/python

import qctoolkit as qtk
import glob

A = qtk.Molecule('data/A.xyz')
B = qtk.Molecule('data/B.xyz')
AB = qtk.Molecule('data/AB.xyz')
p = qtk.Molecule('data/periodic_algaas.cyl')

B.write()
