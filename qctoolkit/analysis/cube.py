# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import qctoolkit as qtk
import re, copy
import qctoolkit.utilities as ut
import read_cube as rq
import write_cube as wq
import esp_point as ESP_c
import esp_cube as ESP_cube_c
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
  dV = x[2] * y[2] * z[2]
  # already in Angstrom?
  shape = (len(x[-1]), len(y[-1]), len(z[-1]))

  zcoord = np.array([[1., 0., 0., 0.]])
  data = c.reshape(*shape, order = 'F') / dV
  data_ref = data.ravel()
  data = data_ref.reshape(*shape)

  return data, zcoord, grid

def read_gaussian(fchk, **kwargs):
  if 'name' in kwargs:
    cube = kwargs['name']
    root, ext = os.path.splitext(cube)
    if ext != '.cube':
      print ext
      qtk.warning("extension .cube is required for many " + \
                  "visualization programs")
  else:
    root, ext = os.path.splitext(fchk)
    cube = root + '.cube'
  qtk.progress("CUBE", "writing file %s\n" % cube)
  if 'flag' not in kwargs:
    flag = 'density=scf'
  else:
    flag = kwargs['flag']
  if 'grid' in kwargs:
    grid = kwargs['grid']
    cmd = '%s 1 %s %s %s -1' % (qtk.gaussian_cubegen_exe, 
                                flag, fchk, cube)
    run = sp.Popen(cmd, shell=True, stdin=sp.PIPE)
    for i in range(len(grid)):
      vec = grid[i]
      if i == 0:
        # for formated output
        msg = '-1 %f %f %f\n' % (vec[1], vec[2], vec[3])
      elif i == 1:
        # for Bohr as grid unit
        msg = '%d %f %f %f\n' % (-vec[0], vec[1], vec[2], vec[3])
      else:
        msg = '%d %f %f %f\n' % (vec[0], vec[1], vec[2], vec[3])
      run.stdin.write(msg)
  else:
    cmd = '%s 1 %s %s %s' % (qtk.gaussian_cubegen_exe, 
                             flag, fchk, cube)
    run = sp.Popen(cmd, shell=True, stdin=sp.PIPE)
  run.stdin.flush()
  run.communicate()
  run.wait()
  q = qtk.CUBE(cube, format='cube')
  zcoord = np.hstack([q.molecule.Z[:, np.newaxis], q.molecule.R])
  zcoord[:,1:4] = zcoord[:,1:4] / 0.529177249
  return q.data, zcoord, q.grid

class CUBE(object):
  """
  read Gaussian CUBE file into numpy 3D array
  together with grid setting and strucutre
  external C module is implemented to access *.cube file
  pointwise addition, substraction, multiplication and division
  are implemented for wavefunction/density analysis
  """
  def __init__(self, cube_file = None, **kwargs):
    if cube_file:
      if not os.path.exists(cube_file):
        ut.exit("CUBE file:%s not found" % cube_file)
      self.path, self.name = os.path.split(cube_file)
      if 'format' not in kwargs:
        ext = os.path.splitext(self.name)[1]
        if ext == '.cube':
          kwargs['format'] = 'cube'
        elif ext == '.casino' or ext == '.dat':
          kwargs['format'] = 'casino'
        elif self.name == 'CHGCAR' or ext =='.vasp':
          kwargs['format'] = 'vasp'
        elif ext == '.fchk':
          kwargs['format'] = 'gaussian'
        else:
          qtk.exit("unknown extension %s" % ext)
      if kwargs['format'] == 'cube':
        self.data, self.zcoord, self.grid, self.coords\
          = rq.read_cube(cube_file)
      elif kwargs['format'] == 'vasp':
        self.data, self.zcoord, self.grid = read_vasp(cube_file)
      elif kwargs['format'] == 'casino':
        self.data, self.zcoord, self.grid = read_casino(cube_file)
      elif kwargs['format'] == 'gaussian':
        self.data, self.zcoord, self.grid = \
          read_gaussian(cube_file, **kwargs)

      self.molecule = qtk.Molecule()
      self.shape = self.data.shape
      self.molecule.R = self.zcoord[:,1:4]*0.529177249
      self.molecule.Z = self.zcoord[:,0]
      self.molecule.N = len(self.zcoord)
      self.molecule.type_list = [ut.Z2n(z) for z in self.molecule.Z]
      self.interp = None
  
      def vec(i):
        return self.grid[i,1:]
    
      self.dV = np.dot(vec(1), np.cross(vec(2), vec(3)))
      self.V = self.dV*self.grid[1,0]*self.grid[2,0]*self.grid[3,0]

    else:
      self.grid = np.zeros([4, 4])
      self.molecule = qtk.Molecule()


  def __repr__(self):
    if np.sum(self.grid) > 0:
      return '\nCUBE mesh:\n' + str(self.grid) +\
             '\n\nmolecule:\n' +\
             str(np.hstack([self.molecule.Z[:, np.newaxis], 
                            self.molecule.R]))
    else:
      return "None\n"

  def __getitem__(self, key):
    new = copy.deepcopy(self)
    new.data = new.data[key]
    corner = new.grid[0, 1:]
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

  def __call__(self, *coord, **kwargs):
    if 'unit' not in kwargs:
      unit = 'angstrom'
    else:
      unit = kwargs['unit'].lower()
      if unit not in ['angstrom', 'bohr']:
        qtk.warning('unit %s not reconized, set to Bohr' % unit)
        unit = 'bohr'

    R = np.atleast_2d(np.array(coord))
    if len(coord) == 3:
      try:
        x, y, z = [float(c) for c in coord]
        R = np.atleast_2d(np.array([x,y,z]))
      except TypeError:
        pass

    if unit == 'angstrom':
      R = [r / 0.529177249 for r in R]

    if not self.interp:
      self.interpolate()

    if len(R) == 1:
      return self.interp(R)[0]
    else:
      return self.interp(R)

  def grid_points(self):
    X, Y, Z = self.meshgrid()
    return np.array(zip(X.ravel(), Y.ravel(), Z.ravel()))

  def build(self, molecule, grid, data):
    self.molecule = molecule
    self.grid = grid
    self.data = data

    def vec(i):
      return self.grid[i,1:]

    self.dV = np.dot(vec(1), np.cross(vec(2), vec(3)))
    self.V = self.dV*self.grid[1,0]*self.grid[2,0]*self.grid[3,0]

  def interpolate(self, **kwargs):

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
    interp = RGI((xs, ys, zs), self.data, method = method, 
                 bounds_error=False, fill_value=0)
    self.interp = interp
    return interp

  def remesh(self, other, **kwargs):
    # maybe reimplemented via scipy.interpolate.rbf
    self = copy.deepcopy(self)

    xo = linepoints(other, 0)
    yo = linepoints(other, 1)
    zo = linepoints(other, 2)
    X, Y, Z = np.meshgrid(xo, yo, zo, indexing='ij')
    points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

    interp = self.interpolate()
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
    cmd = '%s 1 density=scf %s %s -1' % (qtk.gaussian_cubegen_exe, 
                                         fchk, cube)
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
    data = self.data
    grid = self.grid
    N = self.molecule.N
    Z = self.molecule.Z.reshape(N, 1)
    structure = np.hstack([Z, self.molecule.R * 1.889725989])
    if os.path.exists(out):
      ut.prompt("output file:%s exist overwrite?" % out)
    wq.write_cube(out, grid, structure, data)

  def integrate(self, **kwargs):
    if 'power' in kwargs:
      power = kwargs['power']
    else:
      power = 1
    if power > 0:
      return np.sum(np.ravel(self.data**power)) * self.dV
    else:
      return np.sum(np.ravel(abs(self.data))) * self.dV

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
      def getList(i):
        lim = self.grid[0, 1] + (self.grid[1, 0] - 1) * self.grid[1, 1]
        return np.linspace(self.grid[0, i], lim, self.grid[i, 0])
      x = getList(1)
      y = getList(2)
      z = getList(3)
      ut.warning("CUBE.meshgrid returns in the unit of Bohr!")
      return np.meshgrid(x,y,z, indexing='ij')
    else:
      ut.exit("meshgrid only works for cubical system")

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
      if 'no_show' in kwargs and kwargs['no_show']:
        pass
      else:
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
    _grid = copy.deepcopy(self.grid)
    _grid[:, 1:4] = _grid[:, 1:4] * 0.529177249
    O = _grid[0,1:4]
    _max = []
    for i in range(1, 4):
      _max.append(_grid[i,i] * _grid[i,0] + O[i-1])
    _max = np.array(_max)
    _axis = []
    _axis_text = ['$x$', '$y$', '$z$']
    _label = []
    _coord = []
    _range = []
    _origin = []

    for i in range(3):
      line = np.linspace(O[i], _max[i], _grid[i+1, 0])
      if i != axis:
        _axis.append(line)
        _label.append(_axis_text[i])
        _coord.append(list(self.molecule.R[:,i]))
        _range.append([_grid[0, 1+i], _max[i]])
        _origin.append(O[i])
      else:
        if not loc:
          loc = np.argmin(abs(line - level))
        loc_coord = (O[i] + loc * _grid[i+1, i+1])
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
    if 'no_show' in kwargs and kwargs['no_show']:
      pass
    else:
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
      if 'title' in plkwargs:
        title = plkwargs.pop('title')
      else:
        title = None
      fig = plt.figure(name)
      ax = fig.add_subplot(111)
      CS = plt.contour(X, Y, Z, levels, *plargs, **plkwargs)
      CB = plt.colorbar(CS, shrink=0.8, extend='both')
      if title is not None:
        plt.title(title)
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

  def getDipole(self, Z_flag='no_core', component='full'):
    if Z_flag == 'no_core':
      Z = [qtk.ve_list[qtk.Z2n(z)] for z in self.molecule.Z]
    else:
      Z = self.molecule.Z
    if abs(self.integrate() - sum(Z)) > 0.2:
      qtk.setting.quiet = False
      qtk.warning("net charge is not zero! check Z_flag='all'")
    pQ = np.array([sum(Z * self.molecule.R[:,i]) 
      for i in range(3)]) * 1.8897259885789
    qtk.setting.quiet = True
    mg = self.meshgrid()
    qtk.setting.quite = False
    pq = [(self.data * mg[i]).sum()*self.dV for i in range(3)]
    if component == 'full':
      return pQ - pq
    elif component == 'nuc':
      return pQ
    elif component == 'ele':
      return pq

  def range_1D(self):
    qrange =  [self.grid[i,i] * (self.grid[i,0]-1) + self.grid[0,i]
      for i in range(1,4)]
    x, y, z = [np.linspace(
      self.grid[0,i+1], qrange[i], self.grid[i+1,0])
      for i in range(3)]
    return x, y, z

  def ESP(self, coord=None, **kwargs):
    """
    method for electron density
    Note: CUBE file is assumed to be orthorohmbic
    """
    data = self.data
    grid = self.grid
    if 'molecule' not in kwargs:
      mol = self.molecule
    else:
      try:
        mol = copy.deepcopy(qtk.toMolecule(kwargs['molecule']))
      except:
        qtk.exit("error when loading molecule:%s" % \
                 str(kwargs['molecule']))
    N = mol.N

    Q = self.integrate()
    Z_sum = sum(mol.Z)
    ne_diff = abs(Q - Z_sum)
    ve_diff = abs(Q - mol.getValenceElectrons())
    if min(ne_diff, ve_diff) > 1E-2:
      qtk.warning("charge not conserved... ESP is wrong!")
      qtk.warning("charge integrate to %.2f, " % Q + \
                  "while nuclear charge adds to %.2f" % Z_sum)
    if ve_diff < ne_diff:
      Z = [qtk.n2ve(qtk.Z2n(z)) for z in mol.Z]
      Z = np.array(Z).reshape(N, 1)
    else:
      Z = mol.Z.reshape(N, 1)
    structure = np.hstack([Z, mol.R * 1.889725989])
    if coord is not None:
      x, y, z = np.array(coord) * 1.889725989
      V = ESP_c.esp_point(grid, structure, data, x, y, z)
      return V
    else:
      qtk.warning("known issue: unidentifed memory leak")
      out = copy.deepcopy(self)
      out.molecule = mol
      out.data = np.nan_to_num(
        ESP_cube_c.esp_cube(grid, structure, data)
      )
      return out

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

  def __radd__(self, other):
    return self.__add__(other)

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

  def __rsub__(self, other):
    return self.__sub__(other)

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

  def __rmul__(self, other):
    return self.__mul__(other)

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

  def __rdiv__(self, other):
    if isinstance(other, CUBE):
      return self.__div__(other)
    else:
      _out = copy.deepcopy(self)
      _out.data = other / _out.data
      return _out

  def __pow__(self, other):
    if isinstance(other, CUBE):
      ut.exit("power operation is not allowed between CUBE objects")
    else:
      _out = copy.deepcopy(self)
      _out.data = _out.data ** other
      return _out

