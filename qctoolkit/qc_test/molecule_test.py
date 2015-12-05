#!/usr/bin/python

import qctoolkit as qtk
import glob

A = qtk.Molecule( 'data/molecules/A.xyz')
B = qtk.Molecule( 'data/molecules/B.xyz')
AB = qtk.Molecule('data/molecules/AB.xyz')
p = qtk.Molecule('data/molecules/periodic_algaas.cyl')

#print 'yo'
#print B
#B.write_xyz()
#B.align([0,1,0])
#B.write_xyz()
#B.align([0,0,1])
#B.write_xyz()
#
#test = B.principalAxes()
#print test[1][0]

#A.write_xyz()
#A.align([1,1,1])
#A.write_xyz()

A.align()
A.write()
