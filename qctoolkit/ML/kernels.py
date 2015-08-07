import numpy as np

class Gaussian(object):
  def __init__(self, s):
    self._sigma = s
  def evaluate(self, x1, x2):
    distance = np.linalg.norm(np.array(x1)-np.array(x2))
    return np.exp(-0.5*(distance/self._sigma)**2)
    
#def Gaussian(x1, x2, *args, **kwargs):
#  #sigma = kwargs['sigma']
#  sigma = args[0]
#  distance = np.linalg.norm(np.array(x1)-np.array(x2))
#  return np.exp(-0.5*(distance/self._sigma)**2)
