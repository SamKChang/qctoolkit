#!/usr/bin/python

import qctoolkit as qtk
import glob


files = sorted(glob.glob('data/molecules/*'))
mols = []
for i in files:
  print i
  mols.append(qtk.QMInp(i, program='vasp'))#, nuclear_charges=[[1,2,'H.test']]))

#mols[0].setAtom(1,element='C')
#out = mols[3].run(threads=4, save_restart=True, fix_molecule=False)
#print out
mols[3].write()
#out = mols[3].run(threads=4, save_restart=True)
#print out.n_basis
#print out.occupation
#print out.mo_eigenvalues
#print out.mo_coefficients
#print out.nuclear_repulsion
