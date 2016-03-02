import qctoolkit as qtk
import numpy as np
from dlist_2 import dlist_2 as dl2
from dlist_1 import dlist_1 as dl1
from vacf import vacf as vacf_c

class GenericMDInput(object):
  def __init__(self, molecule, **kwargs):
    self.molecule = qtk.toMolecule(molecule)
    if 'temperature' not in kwargs:
      kwargs['temperatur'] = 298
    if 'temperature_tolerance' not in kwargs:
      kwargs['temperature_tolerance'] = 50
    if 'md_sample_period' not in kwargs:
      kwargs['md_sample_period'] = 50
    self.setting = kwargs


class GenericMDOutput(object):
  def __init__(self, out_dir=None, **kwargs):
    self.N = 0
    self.position = None
    self.velocity = None
    self.time = []
    self.type_list = []
    self.cell = None
    self.t_unit = 0
    self.t_step = 0
    self.l_unit = 0

  def __repr__(self):
    return str(self.position)

  def vacf(self, type_list=None, **kwargs):
    """
    velocity autocorrelation function
    """
    if type_list:
      if type(type_list) is not list: type_list = [type_list]
      if type(type_list[0]) is str:
        new = [i for i in range(self.N)\
               if self.type_list[i] in type_list]
        type_list = new
    else: type_list = range(self.N)

    velocity = self.velocity
    size_t, size_n, _ = velocity.shape
    size = size_t * size_n * 3
    flatVel = list(velocity.reshape([size]))
    va, vax, vay, vaz = vacf_c(flatVel, size_t, size_n, type_list)
    if 'components' in kwargs and kwargs['components']:
      return va, vax, vay, vaz
    else:
      return va

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

    if 'dr' not in kwargs:
      kwargs['dr'] = 0.005;

    if np.sum(abs(
      self.cell - np.diag(np.diag(self.cell)))) != 0:
      qtk.exit("radial distribution is only implemented "+\
               "for orthorhmbic cell")

    if 't_start' not in kwargs:
      kwargs['t_start'] = 0

    def distance_list(list1, list2):
      traj = self.position[kwargs['t_start']:]
      size_t, size_n, _ = traj.shape
      size = size_t * size_n * 3
      flatTraj = list(traj.reshape([size]))
      cell = list(np.diag(self.cell))
      if list1 == list2:
        return dl1(flatTraj, size_t, size_n, list1, list2, 
                   cell, kwargs['dr'])
      else:
        return dl2(flatTraj, size_t, size_n, list1, list2, 
                   cell, kwargs['dr'])
        
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
