#!/usr/bin/python

import qctoolkit as qtk
import glob


files = sorted(glob.glob('data/molecules/*'))
mols = []
for i in files:
  mols.append(qtk.QMInp(i, program='cpmd', scf_step=1))

#mols[0].setAtom(1,element='C')
out = mols[0].run()

print out
