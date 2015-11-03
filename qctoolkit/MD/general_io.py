import mdtraj as md
import numpy as np
from dlist_2 import dlist_2 as dl2
from dlist_1 import dlist_1 as dl1

class GenericMDOutput(object):
  def __init__(self, traj, topology):
    self.data = md.load(traj, top=topology)
    self.N = self.data.n_atoms
    self.time = self.data.time
    self.slice = self.data.slice

  def __repr__(self):
    return str(self.data.xyz)

  def gr(self, type1=None, type2=None, **kwargs):
    """
    Radial distribution funciton
    1) 
    Default: calculate all pairwise distances

    2)
    If one atom type specified: calculate only for one type
    e.g: g_OO for water

    3)
    If two atom types specified: calculate distances between
    two specified atom types
    """

    if 't_start' not in kwargs:
      kwargs['t_start'] = 0

    def distance_list(list1, list2):
      traj = self.data.xyz[kwargs['t_start']:]
      size_t, size_n, _ = traj.shape
      size = size_t * size_n * 3
      flatTraj = list(traj.reshape([size]))
      if list1 == list2:
        return dl1(flatTraj, size_t, size_n, list1, list2)
      else:
        return dl2(flatTraj, size_t, size_n, list1, list2)
        
    # the case for two atom types specifed
    if type2:
      # if two specified types differ
      if type1 != type2:
        list1 = [i for i in range(self.N) \
          if self.type_list[i] == type1]
        list2 = [i for i in range(self.N) \
          if self.type_list[i] == type2]
        return distance_list(list1, list2)
      else:
        list1 = [i for i in range(self.N) \
          if self.type_list[i] == type1]
        return distance_list(list1, list1)
    # for the case only one type specified
    elif type1:
      list1 = [i for i in range(self.N) \
        if self.type_list[i] == type1]
      return distance_list(list1, list1)
    # default case:
    else:
      list1 = range(self.N)
      return distance_list(list1, list1)
