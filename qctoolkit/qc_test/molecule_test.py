#!/usr/bin/python

import qctoolkit as qtk
import glob

A = qtk.Molecule( 'data/molecules/A.xyz')
B = qtk.Molecule( 'data/molecules/B.xyz')
AB = qtk.Molecule('data/molecules/AB.xyz')
p = qtk.Molecule('data/molecules/periodic_algaas.cyl')

B.write()
