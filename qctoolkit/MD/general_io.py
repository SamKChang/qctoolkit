import mdtraj as md

class GenericMDOutput(object):
  def __init__(self, traj, topology):
    self.data = md.load(traj, top=topology)
    self.N = self.trajectory.n_atoms

  def __repr__(self):
    return str(self.data.xyz)
