#!/usr/bin/python

import qctoolkit as qtk
import glob


files = sorted(glob.glob('data/molecules/*'))
mols = []
for i in files:
  print i
  mols.append(qtk.QMInp(i, program='nwchem'))

#mols[0].setAtom(1,element='C')
#out = mols[0].run()
#print out
mols[3].write()
out = mols[3].run(save_restart=True)
print out

