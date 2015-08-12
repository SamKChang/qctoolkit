import qctoolkit as qtk
import numpy as np
import random

def mc(target_function, ccs_coord, ccs_span, inp_list, **kwargs):
#  if 'max_step' not in kwargs:
#    qtk.exit("'max_step' is missing")
  if 'T' not in kwargs:
    qtk.exit("'T' is missing")
  T = kwargs['T']
  if 'target' not in kwargs:
    qtk.exit("'target' is missing")
  target = kwargs['target']

  E_list = []
  coord_list = []

  itr = 0
  def EcAppend(_new_E, _new_coord):
    E_list.append(_new_E)
    coord_list.append(_new_coord)
    itr += 1
    if itr > 100:
      itr = 30
      print list(zip(E_list,coord_list))
      E_list = E_list[-30:]
      coord_list = coord_list[-30:]

  def E_average():
    length = min(30, len(E_list))
    return sum(E_list[-length:])/float(length)

  def sample(ccs_new_coord):
    if type(inp_list[-1]) == dict:
      kwgs = inp_list[-1]
      args = inp_list[:-1]
      out = target_function(ccs_new_coord, ccs_span, *args, **kwgs)
    else:
      out = target_function(ccs_new_coord, ccs_span, *inp_list)
    diff = out - target
    boltzmann = np.exp(-abs(diff)/float(T))
    rand = random.uniform(0,1)
    if rand >= boltzmann:
      return diff, ccs_new_coord
    else:
      tmp, new_coord = ccs_span.random()
      sample(new_coord)

  EcAppend(sample(ccs_coord), ccs_coord)
  converge = E_average()
  tmp, new_coord = ccs_span.random()
  while converge > cutoff:
    new_E, new_coord = sample(new_coord)
    EcAppend(new_E, new_coord)
