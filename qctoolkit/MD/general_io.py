import mdtraj as md

class GenericMDOutput(object):
  def __init__(self, traj, topology):
    self.trajectory = md.load(traj, top=topology)
    self.N = self.trajectory.n_atoms

  def __repr__(self):
    return str(self.trajectory.xyz)
