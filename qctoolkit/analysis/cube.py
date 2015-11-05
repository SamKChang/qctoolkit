# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import re, copy
import qctoolkit.molecule as geometry
import qctoolkit.utilities as ut
import read_cube as rq
import qctoolkit.setting as setting
import numpy as np
import matplotlib.pyplot as pl

class CUBE(object):
  """
  read Gaussian CUBE file into numpy 3D array
  together with grid setting and strucutre
  external C module is implemented to access *.cube file
  pointwise addition, substraction, multiplication and division
  are implemented for wavefunction/density analysis
  """
  def __init__(self, cube_file):
    self.data, self.zcoord, self.grid = rq.read_cube(cube_file)
    self.structure = geometry.Molecule()
    if(self.grid[0,0] > 0):
      self.structure.R = self.zcoord[:,1:4]*0.529177249
      self.unit = 'Bohr'
    else:
      self.structure.R = self.zcoord[:,1:4]
      self.unit = 'Angstrom'
    self.structure.Z = self.zcoord[:,0]
    self.structure.N = len(self.zcoord)

  def integrate(self, **kwargs):
    self._O = self.grid[0,1:4]
    self._vx = self.grid[1,1:4]
    self._vy = self.grid[2,1:4]
    self._vz = self.grid[3,1:4]
    self._dV = np.linalg.norm(self._vx - self._O)\
              *np.linalg.norm(self._vy - self._O)\
              *np.linalg.norm(self._vz - self._O)
    if 'power' in kwargs:
      power = kwargs['power']
    else:
      power = 1
    return np.sum(np.ravel(self.data**power)) * self._dV

  def isCubic(self):
    grid_vector = self.grid[1:,1:]
    diff = np.linalg.norm(grid_vector, axis=1)\
           - np.sum(grid_vector, axis=1)
    if np.linalg.norm(diff) > 1E-8:
      ut.exit("only orthorhombic grid is supported")
      return False
    else:
      return True

  def meshgrid(self):
    if self.isCubic():
      xi, yi, zi = (self.grid[0, i] for i in range(1,4))
      dx, dy, dz = (self.grid[i, i] for i in range(1,4))
      xf, yf, zf = (self.grid[i, 0]*self.grid[i, i]\
                    for i in range(1,4))
      # convert cube file Bohr to angstrom
      #xr = np.arange(xi, xf, dx)*0.529177249
      #yr = np.arange(yi, yf, dy)*0.529177249
      #zr = np.arange(zi, zf, dz)*0.529177249
      xr = np.arange(xi, xf, dx)
      yr = np.arange(yi, yf, dy)
      zr = np.arange(zi, zf, dz)
      # NOTE: numpy/matlab x-y switch feature...
      #return np.meshgrid(yr, xr, zr)
      ut.warning("CUBE.meshgrid returns in the unit of Bohr!")
      return np.meshgrid(xr, yr, zr, indexing='ij')
    

  def plot(self, **kwargs):
    if self.isCubic():
      if 'axis' in kwargs:
        plot_axis = kwargs['axis']
      else:
        plot_axis = 0
  
      steps = [self.grid[i, i] for i in range(1,4)\
               if not i == plot_axis+1]
      axes = [i for i in range(3) if not i==plot_axis]
      i = plot_axis+1
      xi, xf, dx = (self.grid[0, i], 
                    self.grid[i, 0]*self.grid[i,i]+self.grid[0,i],
                    self.grid[i, i])
      xout = np.arange(xi, xf, dx)*0.529177249
      yout = np.sum(self.data, axis=tuple(axes))*np.prod(steps)
      if not 'no_show' in kwargs:
        if 'name' in kwargs:
          name = kwargs['name']
        else:
          name = 'density plot'
        if 'plargs' in kwargs:
          plargs = kwargs['plargs']
        else:
          plargs = []
        if 'plkwargs' in kwargs:
          plkwargs = kwargs['plkwargs']
        else:
          plkwargs = {}
  
        pl.figure(name)
        pl.plot(xout, yout, *plargs, **plkwargs)
        if 'legend' in kwargs:
          pl.legend([kwargs['legend']])
      return xout, yout

  def contour(self, axis, level):
    pass

  def shift(self, vector):
    vectorb = np.array(vector) / 0.529177249
    self.grid[0][1:] = self.grid[0][1:] + vectorb
    self.structure.shift(np.array(vector))

  def __add__(self, other):
    if isinstance(other, CUBE):
      _grid = self.grid[1:] - other.grid[1:]
      if(abs(sum(sum(_grid))) < 10**-7 ):
        _out = copy.deepcopy(self)
        _out.data = self.data + other.data
        return _out
      else:
        ut.exit("ERROR from qmout.py->CUBE: " +\
                 "grid mesh does not match")
    else:
      _out = copy.deepcopy(self)
      _out.data = _out.data + other
      return _out

  def __sub__(self, other):
    if isinstance(other, CUBE):
      _grid = self.grid[1:] - other.grid[1:]
      if(abs(sum(sum(_grid))) < 10**-7 ):
        _out = copy.deepcopy(self)
        _out.data = self.data - other.data
        return _out
      else:
        ut.exit("ERROR from qmout.py->CUBE: " +\
                 "grid mesh does not match")
    else:
      _out = copy.deepcopy(self)
      _out.data = _out.data - other
      return _out

  def __mul__(self, other):
    if isinstance(other, CUBE):
      _grid = self.grid[1:] - other.grid[1:]
      if(abs(sum(sum(_grid))) < 10**-7 ):
        _out = copy.deepcopy(self)
        _out.data = np.multiply(self.data, other.data)
        return _out
      else:
        ut.exit("ERROR from qmout.py->CUBE: " +\
                 "grid mesh does not match")
    else:
      _out = copy.deepcopy(self)
      _out.data = _out.data * other
      return _out

  def __div__(self, other):
    if isinstance(other, CUBE):
      _grid = self.grid[1:] - other.grid[1:]
      if(abs(sum(sum(_grid))) < 10**-7 ):
        _out = copy.deepcopy(self)
        _out.data = np.divide(self.data, other.data)
        return _out
      else:
        ut.exit("ERROR from qmout.py->CUBE: " +\
                 "grid mesh does not match")
    else:
      _out = copy.deepcopy(self)
      _out.data = _out.data / other
      return _out
