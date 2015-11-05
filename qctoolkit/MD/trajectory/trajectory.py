import mdtraj as md

class Trajectory(object):
  def __init__(self, traj=None, topology=None):
    if traj and topology:
      self.data = md.load(traj, top=topology)
      self.N = self.data.n_atoms
      self.time = self.data.time
      self.slice = self.data.slice
    else:
      self.data = None
      self.N = 0
      self.time = []
    self.cell = []

