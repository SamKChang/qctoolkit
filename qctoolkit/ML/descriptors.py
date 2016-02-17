import qctoolkit as qtk
import numpy as np
import coulomb_matrix as cm
import re

class CoulombMatrix(object):
  def __init__(self, *args, **kwargs):
    if len(args) == 1 and args[0] > 0:
      self.base_dim = args[0]
    elif len(args) == 0 or args[0] == 0:
      fxyz = open(kwargs['xyz'], "r")
      self.base_dim = int(fxyz.readline())
      fxyz.close()
    
    if 'xyz' in kwargs:
      self.data = cm.coulomb_matrix(kwargs['xyz'], 
                                         self.base_dim)
      self.name = re.sub('.*\/','',kwargs['xyz'])
      self.name = re.sub('\.xyz','',self.name)
    else:
      qtk.exit("CoulombMatrix: input mode is not specified")

    # set row norm matrix
    NORM = np.vstack([sum(self.data), range(self.base_dim)]).T
    # sort NORM by row norm
    NORM=NORM[NORM[:,0].argsort()]
    # reverse array order
    NORM=NORM[::-1]
    # extract new row order
    sortRow=NORM[:,1].astype(int)
    # rearrange Coulomb matrix
    self.data = self.data[:,sortRow][sortRow]

  
  def getVector(self):
    return self.data.reshape(self.base_dim**2), self.name
  
