#!/usr/bin/python

import qctoolkit as qtk
import glob


files = sorted(glob.glob('data/molecules/*'))
mols = []
for i in files:
  mols.append(qtk.QMInp(i, program='vasp'))

#mols[0].setAtom(1,element='C')
mols[0].write()
#mols[0].run()
