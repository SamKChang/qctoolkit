# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import qctoolkit as qtk
import re, copy
import qctoolkit.molecule as geometry
import qctoolkit.utilities as ut
import read_cube as rq
import write_cube as wq
import qctoolkit.setting as setting
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import subprocess as sp
mpl.use('Agg')

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
    self.path, self.name = os.path.split(cube_file)
    self.data, self.zcoord, self.grid, self.coords\
      = rq.read_cube(cube_file)
    self.coords = self.coords * 0.529177249
    self.molecule = geometry.Molecule()
    self.shape = self.data.shape
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

  def __getitem__(self, key):
    new = copy.deepcopy(self)
    new.data = new.data[key]
    corner = new.grid[0, 1:]
    print key
    for i in range(len(key)):
      k = key[i]
      v = self.grid[i+1, 1:]
      if type(k) is int:
        start = k
        end = 1
      else:
        start = k.start
        end = k.stop
      if start is None:
        start = 0
      if end is None:
        end = int(self.grid[1+i, 0])
      corner = corner + v*start
      new.grid[1+i, 0] = end
    new.grid[0, 1:] = corner
    return new

  def asGaussianTemplate(self, gcube, **kwargs):
    cube_name_list = gcube.name.split('.')
    fchk_name = '.'.join(cube_name_list[:-1]) + '.fchk'
    fchk = os.path.abspath(os.path.join(gcube.path, fchk_name))
    path, _ = os.path.split(fchk)
    if not os.path.exists(fchk):
      ut.exit("gaussian fchk file:%s not found" % fchk)
    cube = 'gcube_tmp_' + str(id(self)) + '.cube'
    cube = os.path.join(path, cube)
    cmd = '%s 1 density=scf %s %s -1' % (qtk.cubegen_exe, fchk, cube)
    run = sp.Popen(cmd, shell=True, stdin=sp.PIPE)
    for i in range(len(self.grid)):
      vec = self.grid[i]
      if i == 0:
        # for formated output
        msg = '-1 %f %f %f\n' % (vec[1], vec[2], vec[3])
      elif i == 1:
        # for Bohr as grid unit
        msg = '%d %f %f %f\n' % (-vec[0], vec[1], vec[2], vec[3])
      else:
        msg = '%d %f %f %f\n' % (vec[0], vec[1], vec[2], vec[3])
      run.stdin.write(msg)
    run.stdin.flush()
    run.communicate()
    run.wait()
    return qtk.CUBE(cube)

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
  
        plt.figure(name)
        plt.plot(xout, yout, *plargs, **plkwargs)
        if 'legend' in kwargs:
          plt.legend([kwargs['legend']])
      return xout, yout

  def filter(self, cutoff=0.001):
    m = np.max(self.data)
    self.data[np.abs(self.data) < m*cutoff] = 0

  def contour(self, axis=0, **kwargs):
    if 'slice' not in kwargs:
      loc = None
    else:
      loc = kwargs['slice']
    if not loc:
      level = self.molecule.getCenterOfMass()[axis]
      _text = ['x', 'y', 'z']
      ut.report("CUBE", 
                "center of mass on %s-axis:" % _text[axis], level)
      level = level / 0.529177249
    O = self.grid[0,1:4]
    _max = []
    for i in range(1, 4):
      _max.append(self.grid[i,i] * self.grid[i,0])
    _max = np.array(_max)
    _axis = []
    _axis_text = ['$x$', '$y$', '$z$']
    _label = []
    for i in range(3):
      line = np.linspace(O[i], _max[i], self.grid[i+1, 0])
      if i != axis:
        _axis.append(line * 0.529177249)
        _label.append(_axis_text[i])
      else:
        if not loc:
          loc = np.argmin(abs(line - level))
    if axis == 0:
      try:
        Z = self.data[loc, :, :]
      except:
        qtk.exit("failed when accessing cube data")
    elif axis == 1:
      try:
        Z = self.data[:, loc, :]
      except:
        qtk.exit("failed when accessing cube data")
    elif axis == 2:
      try:
        Z = self.data[:, :, loc]
      except:
        qtk.exit("failed when accessing cube data")
        
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
    plt.figure(name)
    CS = plt.contour(X, Y, Z, 10, *plargs, **plkwargs)
    CB = plt.colorbar(CS, shrink=0.8, extend='both')
    #plt.clabel(CS, fontsize=9, inline=1)
    plt.xlabel(_label[0] + r" [$\rm \AA$]", fontsize=15)
    plt.ylabel(_label[1] + r" [$\rm \AA$]", fontsize=15)
    plt.axes().set_aspect('equal', 'datalim')
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
