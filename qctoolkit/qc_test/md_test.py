#!/usr/bin/python

import qctoolkit as qtk
import qctoolkit.MD as qmd
import pylab as p
import numpy as np
import matplotlib.pyplot as plt

out = 'data/mdout/cpmd/'

traj = qmd.MDOut(out)
#print traj.position
#print traj.velocity
#print traj.cell

out = traj.vacf()
#for x0,y0 in zip(x, out):
#  print x0,y0
#
#print 'pos(11): ',
#print traj.position[0,11,:]
#print 'pos(22): ',
#print traj.position[0,22,:]

#traj = qmd.xyzOutput('data/mdout/TRAJEC-short.xyz')
#
#out = traj.gr('O')
plt.plot(out)
#n, bins, patches = p.hist(out, bins=100)
#np.savetxt('traj.dat', out)
p.show()
