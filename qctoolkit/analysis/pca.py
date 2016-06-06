from __future__ import division
import qctoolkit as qtk
import numpy as np
import copy
dot = np.dot

class PCA:
  def __init__( self, A, **kwargs):
    if 'fraction' not in kwargs:
      fraction = 0.9
    else:
      fraction = kwargs['fraction']
    if 'axis' not in kwargs:
      axis = 0
    else:
      axis = kwargs['axis']
    A = copy.deepcopy(A)
    self.std = A.std(axis=axis)
    self.mean = A.mean(axis=axis)
    if 'scale' in kwargs and kwargs['scale']:
      qtk.report('PCA', 'rescaling data')
      A /= self.std
    qtk.report("PCA", "centering data")
    A -= self.mean

    assert 0 <= fraction <= 1
        # A = U . diag(d) . Vt, O( m n^2 ), lapack_lite --
    self.U, self.d, self.Vt = np.linalg.svd( A, full_matrices=False )
    assert np.all( self.d[:-1] >= self.d[1:] )  # sorted
    self.eigen = self.d**2
    self.sumvariance = np.cumsum(self.eigen)
    self.sumvariance /= self.sumvariance[-1]
    self.npc = np.searchsorted( self.sumvariance, fraction ) + 1
    self.dinv = np.array([ 1/d if d > self.d[0] * 1e-6  else 0
                            for d in self.d ])
    self.data = A

  def pc( self ):
    """ e.g. 1000 x 2 U[:, :npc] * d[:npc], to plot etc. """
    n = self.npc
    return self.U[:, :n] * self.d[:n]

  # These 1-line methods may not be worth the bother;
  # then use U d Vt directly --

  def vars_pc( self, x ):
    n = self.npc
    # 20 vars -> 2 principal
    return self.d[:n] * dot( self.Vt[:n], x.T ).T  

  def pc_vars( self, p ):
    n = self.npc
    # 2 PC -> 20 vars
    return dot( self.Vt[:n].T, (self.dinv[:n] * p).T ) .T  

  def pc_obs( self, p ):
    n = self.npc
    # 2 principal -> 1000 obs
    return dot( self.U[:, :n], p.T )  

  def obs_pc( self, obs ):
    n = self.npc
    # 1000 obs -> 2 principal
    return dot( self.U[:, :n].T, obs ) .T  

  def obs( self, x ):
    # 20 vars -> 2 principal -> 1000 obs
    return self.pc_obs( self.vars_pc(x) )  

  def vars( self, obs ):
    # 1000 obs -> 2 principal -> 20 vars
    return self.pc_vars( self.obs_pc(obs) )  
