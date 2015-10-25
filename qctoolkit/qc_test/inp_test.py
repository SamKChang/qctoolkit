#!/usr/bin/python

import qctoolkit as qtk
import glob


files = sorted(glob.glob('data/*'))
mols = []
for i in files:
  print files
  mols.append(qtk.QMInp(i, program='vasp'))

#mols[0].setAtom([2,4,6],string='test.psp')
#mols[0].setAtom([2,3,5],string='test2.psp')
#mols[0].write()

mols[0].write()
#mols[0].run('test')
