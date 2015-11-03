import qctoolkit as qtk
import mdtraj as md
import numpy as np

def MDResult(traj, **kwargs):

  mol = qtk.Molecule(traj)
  mol.write('qtk_tmp_traj.pdb', format='pdb')


#traj1 = md.load('TRAJEC.xyz', top='test.pdb')
#print traj1
traj = md.load('TRAJEC.xyz', top='TRAJEC.pdb')
print traj
