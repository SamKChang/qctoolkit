# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, Gaussian, NwChem, espresso, gamess

import re, sys
import geometry
import numpy as np
import utilities as ut
#import setting
from qctoolkit.io_format import *

class QMInp(setting.QMSetting):
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

  def setChargeMultiplicity(self, charge, multiplicity):
    self.inp.structure.charge = charge
    self.inp.structure.setMultiplicity(multiplicity)

  def setTheory(self, theory):
    self.inp.setting.theory = theory

  def setSCFStep(self, step):
    self.inp.setting.maxstep = step
    self.inp.set_step = True

  def restart(self):
    self.inp.restart = True

  def debug(self):
    self.inp.debug = True

  def write(self, name):
    self.inp.write(name)

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
