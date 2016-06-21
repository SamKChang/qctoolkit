# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import qctoolkit as qtk
import re, copy
import qctoolkit.utilities as ut
import read_cube as rq
import write_cube as wq
import qctoolkit.setting as setting
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import subprocess as sp
from scipy.interpolate import RegularGridInterpolator as RGI

def read_vasp(chg_file):
  with open(chg_file) as cfile:
    content = cfile.readlines()
  grid_1d = np.fromstring(''.join(content[2:5]), dtype=float, sep=' ')
  grid = grid_1d.reshape([3,3])
  type_list = filter(None, content[5].split(' '))[:-1]
  n_list = np.fromstring(content[6], dtype=int, sep=' ')
  Z = []
  for a in range(len(type_list)):
    for n in range(n_list[a]):
      Z.append(qtk.n2Z(type_list[a]))
  Z = np.array(Z)[:, None]
  N = sum(n_list)
  R_scale = np.fromstring(''.join(content[8:8+N]), 
                          dtype=float, sep=' ').reshape([N, 3])
  R = []
  for i in range(N):
    r = np.array([0,0,0])
    for j in range(3):
      r = r + R_scale[i, j] * grid.T[j] / 0.529177249
    R.append(r)
  R = np.array(R)
  zcoord = np.hstack([Z, R])

  V = np.dot(np.cross(grid[0], grid[1]), grid[2]) / (0.529177249)**3
  step = np.fromstring(content[9+N], dtype=int, sep=' ')
  for i in range(3):
    grid[i] = grid[i] / (step[i] * 0.529177249)

  # numpy array.reshape does not change memory order
  # use ravel to unwine new order
  data = np.fromstring(''.join(content[10+N:]), dtype=float, sep=' ')
  data = data.reshape(step, order='F') / V
  data_rav = data.ravel()
  data = data_rav.reshape(step)

  step = step[:, None]
  grid = np.hstack([step, grid])
  grid = np.vstack([[N, 0, 0, 0], grid])

  return data, zcoord, grid

def read_casino(chg_file):
  data = np.loadtxt(chg_file)

  def getRange(np_vector):
    crd = sorted(list(set(list(np_vector))))
    # min, max, step
    return crd[0], crd[-1], crd[1] - crd[0], crd

  x = getRange(data[:,0])
  y = getRange(data[:,1])
  z = getRange(data[:,2])
  c = data[:,3]

  grid = np.array(
    [
      [1, x[0], y[0], z[0]],
      [len(x[-1]), x[2], 0., 0.],
      [len(y[-1]), 0., y[2], 0.],
      [len(z[-1]), 0., 0., z[2]],
    ]
  )
  # already in Angstrom?
  V = (x[1] - x[0]) * (y[1] - y[0]) * (z[1] - z[0])
  shape = (len(x[-1]), len(y[-1]), len(z[-1]))

  zcoord = np.array([[1., 0., 0., 0.]])
  data = c.reshape(*shape, order = 'F') / V
  data_ref = data.ravel()
  data = data_ref.reshape(*shape)

  return data, zcoord, grid

class CUBE(object):
  """
  read Gaussian CUBE file into numpy 3D array
  together with grid setting and strucutre
  external C module is implemented to access *.cube file
  pointwise addition, substraction, multiplication and division
  are implemented for wavefunction/density analysis
  """
  def __init__(self, cube_file, **kwargs):
    if not os.path.exists(cube_file):
      ut.exit("CUBE file:%s not found" % cube_file)
    self.path, self.name = os.path.split(cube_file)
    if 'format' not in kwargs:
      ext = os.path.splitext(self.name)[1]
      if ext == '.cube':
        kwargs['format'] = 'cube'
      elif ext == '.casino':
        kwargs['format'] = 'casino'
      elif self.name == 'CHGCAR' or ext =='.vasp':
        kwargs['format'] = 'vasp'
    if kwargs['format'] == 'cube':
      self.data, self.zcoord, self.grid, self.coords\
        = rq.read_cube(cube_file)
    elif kwargs['format'] == 'vasp':
      self.data, self.zcoord, self.grid = read_vasp(cube_file)
    elif kwargs['format'] == 'casino':
      self.data, self.zcoord, self.grid = read_casino(cube_file)
    self.molecule = qtk.Molecule()
    self.shape = self.data.shape
    self.molecule.R = self.zcoord[:,1:4]*0.529177249
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

  def remesh(self, other, **kwargs):
    # maybe reimplemented via scipy.interpolate.rbf
    self = copy.deepcopy(self)

    def linepoints(cube, i):
      return np.linspace(cube.grid[0, i+1],
                         cube.grid[i+1,0]*cube.grid[i+1,i+1],
                         cube.grid[i+1,0])

    if 'method' in kwargs:
      method = kwargs['method']
    else:
      method = 'linear'

    xs = linepoints(self, 0)
    ys = linepoints(self, 1)
    zs = linepoints(self, 2)
    xo = linepoints(other, 0)
    yo = linepoints(other, 1)
    zo = linepoints(other, 2)
    X, Y, Z = np.meshgrid(xo, yo, zo, indexing='ij')
    points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
    interp = RGI((xs, ys, zs), self.data, method = method, 
                 bounds_error=False, fill_value=0)
    self.data = interp(points)
    self.grid = copy.deepcopy(other.grid)

    return self

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
    out = copy.deepcopy(self)
    out.data[np.abs(self.data) < m*cutoff] = 0
    return out

  def locate(self, coord):
    O = self.grid[0,1:4]
    _max = []
    for i in range(1, 4):
      _max.append(self.grid[i,i] * self.grid[i,0])
    _max = np.array(_max)
    out = []
    for i in range(3):
      line = np.linspace(O[i], _max[i], self.grid[i+1, 0])
      line = line * 0.529177249
      out.append(np.argmin(abs(line - coord[i])))
    return out

  def contour(self, axis=0, **kwargs):
    if 'levels' not in kwargs:
      levels = 15
    else:
      levels = kwargs['levels']
    if 'filter' in kwargs:
      self = self.filter(kwargs['filter'])
    if 'slice' not in kwargs:
      loc = None
    else:
      loc = kwargs['slice']
    if not loc:
      level = self.molecule.getCenterOfMass()[axis]
      _text = ['x', 'y', 'z']
      ut.report("CUBE", 
                "center of mass on %s-axis:" % _text[axis], level)
      level = level
    O = self.grid[0,1:4]
    _max = []
    for i in range(1, 4):
      _max.append(self.grid[i,i] * self.grid[i,0])
    _max = np.array(_max)
    _axis = []
    _axis_text = ['$x$', '$y$', '$z$']
    _label = []
    _coord = []
    _range = []

    for i in range(3):
      line = np.linspace(O[i], _max[i], self.grid[i+1, 0])
      line = line * 0.529177249
      if i != axis:
        _axis.append(line)
        _label.append(_axis_text[i])
        _coord.append(list(self.molecule.R[:,i]))
        _range.append([self.grid[0, 1+i]* 0.529177249, 
                       _max[i]* 0.529177249])
      else:
        if not loc:
          loc = np.argmin(abs(line - level))
        loc_coord = (O[i] + loc * self.grid[i+1, i+1]) * 0.529177249
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
      name = self.name
    if 'plargs' in kwargs:
      plargs = kwargs['plargs']
    else:
      plargs = []
    plkwargs = {}
    if 'plkwargs' in kwargs:
      plkwargs.update(kwargs['plkwargs'])
    fig = plt.figure(name)
    ax = fig.add_subplot(111)
    CS = plt.contour(X, Y, Z, levels, *plargs, **plkwargs)
    CB = plt.colorbar(CS, shrink=0.8, extend='both')
    plt.xlabel(_label[0] + r" [$\rm \AA$]", fontsize=15)
    plt.ylabel(_label[1] + r" [$\rm \AA$]", fontsize=15)
    plt.axes().set_aspect('equal')

    x_list = []
    y_list = []

    def plotElement(i):
      to_plot = True
      for j in range(2):
        if _coord[j][i] > _range[j][1] or _coord[j][i] < _range[j][0]:
          to_plot = False
      if to_plot:
        symbol = self.molecule.type_list[i]
        x = _coord[0][i]
        y = _coord[1][i]
        x_list.append(x)
        y_list.append(y)
        ax.annotate(symbol, xytext=(x+0.02, y+0.02), xy=(0, 0))

    for i in range(self.molecule.N):
      plotElement(i)

    plt.plot(x_list, y_list, ls='', marker='o', color='k')
    x_min, x_max = _range[0]
    y_min, y_max = _range[1]
    plt.xlim(_range[0])
    plt.ylim(_range[1])

    ut.report('CUBE', 'axis:%d, slice:%f' % (axis, loc))
    ut.report("CUBE", "slice coordinate: %f" % loc_coord)
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
