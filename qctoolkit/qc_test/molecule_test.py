#!/usr/bin/python

import qctoolkit as qtk
import glob


files = sorted(glob.glob('data/*'))
mols = []
for i in files:
  print files
  mols.append(qtk.Molecule(i))

mols[0].setAtom(1,element='C')
mols[0].write_xyz()
