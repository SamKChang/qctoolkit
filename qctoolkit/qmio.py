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

class CUBE(object):
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

  def integrate(self, flag):
    self._O = self.grid[0,1:4]
    self._vx = self.grid[1,1:4]
    self._vy = self.grid[2,1:4]
    self._vz = self.grid[3,1:4]
    self._dV = np.linalg.norm(self._vx - self._O)\
              *np.linalg.norm(self._vy - self._O)\
              *np.linalg.norm(self._vz - self._O)
    if flag==1:
      return np.sum(np.ravel(self.data)) * self._dV
    elif flag==2:
      return np.sum(np.ravel(self.data**2)) * self._dV

  def plot(self, axis):
    pass

  def contour(self, axis, level):
    pass

  def __add__(self, other):
    if isinstance(other, CUBE):
      _grid = self.grid - other.grid
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
      _grid = self.grid - other.grid
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
      _grid = self.grid - other.grid
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
      _grid = self.grid - other.grid
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
  def __init__(self, structure_inp, program, info):
    self.program = program
    self.info = info
    self.atom_count = 0

    # take inpur 'program' to choose corresponding format
    if re.match('cpmd', self.program):
      self.inp = cpmd.inp(structure_inp, self.info)

  def setAtom(self, atom_list, atom_string):
    for I in atom_list:
      i = I-1
      self.inp.structure.type_list[i] = atom_string
      self.inp.structure.Z[i] = self.atom_count
    self.inp.atom_list[str(self.atom_count)] = atom_string
    self.atom_count =- 1

  def setCorner(self, corner_coord):
    self.inp.structure.center(-np.array(corner_coord))
    self.inp.set_center = True

  def setCenter(self, center_coord):
    self.inp.setting.center = center_coord
    self.inp.set_center = True

  def setCelldm(self, celldm):
    self.inp.setting.celldm = celldm
    self.inp.set_celldm = True

  def setMargin(self, margin):
    self.inp.setting.margin = margin
    self.inp.set_margin = True

  def setMode(self, mode):
    self.inp.setting.mode = mode
    self.inp.set_mode = True

  def setChargeMultiplicity(self, charge, multiplicity, **kargs):
    self.inp.structure.charge = charge
    self.inp.structure.setMultiplicity(multiplicity, **kargs)

  def setTheory(self, theory):
    self.inp.setting.theory = theory

  def setSCFStep(self, step):
    self.inp.setting.maxstep = step
    self.inp.set_step = True

  def setInitRandom(self):
    self.inp.set_init_random = True

  def restart(self):
    self.inp.restart = True

  def debug(self):
    self.inp.debug = True

  def write(self, name):
    self.inp.write(name)

  def periodic(self):
    self.inp.setting.isolated = False
    self.inp.set_center = False

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

  def getEt(self, name):
    out = open(name, 'r')
