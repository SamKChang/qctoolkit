#!/usr/bin/python

import qctoolkit as qtk
import qctoolkit.MD as qmd

traj = qmd.xyzOutput('data/mdout/TRAJEC.xyz')

print traj
print traj.type_list

traj.gr('H', 'O')


