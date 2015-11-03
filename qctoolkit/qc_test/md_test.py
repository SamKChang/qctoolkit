#!/usr/bin/python

import qctoolkit as qtk
import qctoolkit.MD as qmd
import pylab as p

traj = qmd.xyzOutput('data/mdout/TRAJEC.xyz')

out = traj.gr()
n, bins, patches = p.hist(out)

p.show()
