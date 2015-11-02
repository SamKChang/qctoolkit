import mdtraj as md

class GenericMDOutput(object):
  def __init__(self, traj, topology):
    self.data = md.load(traj, top=topology)
    self.N = self.data.n_atoms

  def __repr__(self):
    return str(self.data.xyz)

  def gr(self, type1=None, type2=None, **kwargs):
    if type2:
      if type1 != type2:
        list1 = [i for i in range(self.N) \
          if self.type_list[i] == type1]
        list2 = [i for i in range(self.N) \
          if self.type_list[i] == type2]
        for i in list1:
          for j in list2:
            print i,j
