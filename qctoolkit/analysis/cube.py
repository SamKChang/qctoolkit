# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import re, copy
import qctoolkit.molecule as geometry
import qctoolkit.utilities as ut
import read_cube as rq
import write_cube as wq
import qctoolkit.setting as setting
import numpy as np
import matplotlib.pyplot as pl
import os

class CUBE(object):
  """
  read Gaussian CUBE file into numpy 3D array
  together with grid setting and strucutre
  external C module is implemented to access *.cube file
  pointwise addition, substraction, multiplication and division
  are implemented for wavefunction/density analysis
  """
  def __init__(self, cube_file):
    if not os.path.exists(cube_file):
      ut.exit("CUBE file:%s not found" % cube_file)
    self.data, self.zcoord, self.grid, self.coords\
      = rq.read_cube(cube_file)
    self.coords = self.coords * 0.529177249
    self.molecule = geometry.Molecule()
    if(self.grid[0,0] > 0):
      self.molecule.R = self.zcoord[:,1:4]*0.529177249
      self.unit = 'Bohr'
    else:
      self.molecule.R = self.zcoord[:,1:4]
      self.unit = 'Angstrom'
    self.molecule.Z = self.zcoord[:,0]
    self.molecule.N = len(self.zcoord)
    self.molecule.type_list = [ut.Z2n(z) for z in self.molecule.Z]

  def __repr__(self):
    return '\nCUBE mesh:\n' + str(self.grid) +\
           '\n\nmolecule:\n' +\
           str(np.hstack([self.molecule.Z[:, np.newaxis], 
                          self.molecule.R]))

  def write(self, out):
    x, y, z = self.data.shape
    data = self.data
    grid = self.grid
    N = self.molecule.N
    Z = self.molecule.Z.reshape(N, 1)
    structure = np.hstack([Z, self.molecule.R * 1.889725989])
    if os.path.exists(out):
      ut.prompt("output file:%s exist overwrite?" % out)
    wq.write_cube(out, grid, structure, data)

  def integrate(self, **kwargs):
    O = self.grid[0,1:4]
    vx = self.grid[1,1:4]
    vy = self.grid[2,1:4]
    vz = self.grid[3,1:4]
    dV = np.linalg.norm(vx - O)\
              *np.linalg.norm(vy - O)\
              *np.linalg.norm(vz - O)
    if 'power' in kwargs:
      power = kwargs['power']
    else:
      power = 1
    if power > 0:
      return np.sum(np.ravel(self.data**power)) * dV
    else:
      return np.sum(np.ravel(abs(self.data))) * dV

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

  def contour(self, axis=0, level=None, **kwargs):
    if not level:
      level = self.molecule.getCenterOfMass()[axis]
    level = level / 0.529177249
    O = self.grid[0,1:4]
    _max = []
    for i in range(1, 4):
      _max.append(self.grid[i,i] * self.grid[i,0])
    _max = np.array(_max)
    _axis = []
    for i in range(3):
      line = np.linspace(O[i], _max[i], self.grid[i+1, 0])
      if i != axis:
        _axis.append(line)
      else:
        loc = np.argmin(abs(line - level))
    if axis == 0:
      Z = self.data[loc, :, :]
    elif axis == 1:
      Z = self.data[:, loc, :]
    elif axis == 2:
      Z = self.data[:, :, loc]
        
    X, Y = np.meshgrid(*tuple(_axis), indexing='ij')
    if 'name' in kwargs:
      name = kwargs['name']
    else:
      name = 'density contour plot'
    if 'plargs' in kwargs:
      plargs = kwargs['plargs']
    else:
      plargs = []
    plkwargs = {}
    if 'plkwargs' in kwargs:
      plkwargs.update(kwargs['plkwargs'])
    pl.figure(name)
    CS = pl.contour(X, Y, Z, 10, *plargs, **plkwargs)
    pl.clabel(CS, fontsize=9, inline=1)
    ut.report('CUBE', 'axis:%d, slice:%f' % (axis, loc))
    return [X, Y, Z]
    

  def shift(self, vector):
    vectorb = np.array(vector) / 0.529177249
    self.grid[0][1:] = self.grid[0][1:] + vectorb
    self.molecule.shift(np.array(vector))

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
