# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import re, sys, copy
import geometry
import numpy as np
import utilities as ut
from qctoolkit.io_format import *
import qctoolkit.read_cube as rq
import pylab as pl

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
        sys.exit("ERROR from qmio.py->CUBE: " +\
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
        sys.exit("ERROR from qmio.py->CUBE: " +\
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
        sys.exit("ERROR from qmio.py->CUBE: " +\
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
        sys.exit("ERROR from qmio.py->CUBE: " +\
                 "grid mesh does not match")
    else:
      _out = copy.deepcopy(self)
      _out.data = _out.data / other
      return _out

class QMInp(object):
  def __init__(self, structure_inp, program, **kwargs):
    self.program = program
    self.atom_count = 0
    if 'info' in kwargs:
      self.info = kwargs['info']
    else: 
      self.info = structure_inp
#    if 'set_charge' in kwargs and kwargs['set_charge']:
#      self.set_charge = True
#    else:
#      self.set_charge = False

    # take input 'program' to choose corresponding format
    if re.match('cpmd', self.program):
      self.inp = cpmd.inp(structure_inp, self.info)
#                          set_charge=self.set_charge)

  def setAtom(self, atom_list, atom_string):
    _tmp = self.atom_count
    if type(atom_list) == int:
      atom_list = [atom_list]
    for I in atom_list:
      i = I-1
      self.inp.structure.type_list[i] = atom_string
      self.inp.structure.Z[i] = self.atom_count
    self.inp.atom_list[str(self.atom_count)] = atom_string
    self.atom_count = _tmp - 1

  def addAtom(self, pp_string, coordinate):
    structure = copy.deepcopy(self.inp.structure)
    structure.R = np.vstack([structure.R, coordinate])
    structure.Z = np.hstack([structure.Z, 1])
    structure.type_list.append('H')
    structure.N += 1
    self.setStructure(structure, reset=False)
    self.setAtom(structure.N, pp_string)

  def setInfo(self, info):
    self.inp.info = info

  def setStructure(self, structure, **kwargs):
    if 'reset' in kwargs:
      reset=kwargs['reset']
    else:
      reset=True
    if reset:
      self.atom_count = 0
    self.inp.structure = copy.deepcopy(structure)

  def setConvergence(self, convergence):
    self.inp.setting.convergence = convergence
    self.inp.setting.set_convergence = True

  def setCorner(self, corner_coord):
    self.inp.setting.center(-np.array(corner_coord))
    self.inp.setting.set_center = True

  def saveDensity(self):
    self.inp.setting.save_density = True

  def setShift(self, shift_coord):
    self.inp.setting.shift = np.array(shift_coord)
    self.inp.setting.set_shift = True
    self.inp.setting.set_margin = True
    self.inp.setting.set_center = False

  def setCenter(self, center_coord):
    self.inp.setting.center = center_coord
    self.inp.setting.set_center = True

  def setCutoff(self, cutoff):
    self.inp.setting.cutoff = cutoff

  def setCelldm(self, celldm):
    self.inp.setting.celldm = celldm
    self.inp.setting.set_celldm = True

  def setKmesh(self, kmesh):
    self.inp.setting.kmesh = kmesh
    self.inp.kpoints = True

  def setMargin(self, margin):
    self.inp.setting.margin = margin
    self.inp.setting.set_margin = True

  def setMode(self, mode):
    self.inp.setting.mode = mode
    self.inp.setting.set_mode = True

  def setChargeMultiplicity(self, charge, multiplicity, **kargs):
    self.inp.setting.charge = charge
    self.inp.setting.multiplicity = multiplicity
    self.inp.structure.setChargeMultiplicity(charge,
                                             multiplicity,
                                             **kargs)
    self.inp.setting.set_charge = True
    self.inp.setting.set_multiplicity = True

  def setTheory(self, theory):
    self.inp.setting.theory = theory

  def setSCFStep(self, step):
    self.inp.setting.maxstep = step
    self.inp.setting.set_step = True

  def setInitRandom(self):
    self.inp.setting.set_init_random = True

  def restart(self):
    self.inp.setting.restart = True

  def debug(self):
    self.inp.setting.debug = True

  def periodic(self):
    self.inp.setting.isolated = False
    self.inp.setting.set_center = False

  def removeAtom(self, index):
    self.inp.structure.remove_atom(index)

  def isolateAtoms(self, indices):
    self.inp.structure.isolate_atoms(indices)

  def write(self, *args, **kwargs):
    mul = self.inp.structure.multiplicity
    chg = self.inp.structure.charge
    ve = np.vectorize(ut.n2ve)
    nve = sum(ve(self.inp.structure.type_list)) - chg
    #if mul % 2 != (np.sum(self.inp.structure.Z) + chg) % 2:
    if mul % 2 != nve % 2:
      self.inp.write(*args, **kwargs)
    else:
      msg = "Multiplicity %d " % mul + \
            "and %d valence electrons " % nve +\
            "\n(with charge %3.1f) " % float(chg) +\
            "are not compatible"
      ut.exit(msg)

class QMOut(object):
  def __init__(self, qmout, program):
    #self.name = re.sub("\..*", "", qmout)
    #self.file_name = qmout
    self.program = program
    if (re.match('cpmd', program)):
      self.out = cpmd.out(qmout)
    self.Ehartree = self.out.Et
    self.info = self.out.info
    self.Et = self.Ehartree
    self.SCFStep = self.out.SCFStep

  def Ha2ev(self):
    self.Et = self.Ehartree * 27.211396132

  def Ha2kcal(self):
    self.Et = self.Ehartree * 627.509469

  #def getEt(self, name):
  #  out = open(name, 'r')
