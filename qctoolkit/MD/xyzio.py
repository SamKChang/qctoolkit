import qctoolkit as qtk
import re, os
from general_io import GenericMDOutput

class xyzOutput(GenericMDOutput):
  def __init__(self, traj):
    self.molecule = qtk.Molecule(traj)
    stem, ext = os.path.splitext(traj)
    pdb = re.sub('\.xyz', '', traj) + '_qtk_tmp.pdb'
    self.molecule.write(pdb, format='pdb')
    GenericMDOutput.__init__(self, traj, pdb)
    os.remove(pdb)

    self.type_list = self.molecule.type_list
    self.Z = self.molecule.Z
