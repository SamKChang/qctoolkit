# interface wrapper of I/O format and settings
# for general QM calculations
# QM programs include: 
#  CPMD, VASP, Gaussian, NwChem, espresso, gamess

import re, copy, os
import qctoolkit.geometry as geometry
import numpy as np
import qctoolkit.utilities as ut
import cpmd
import vasp
import qctoolkit.read_cube as rq
import qctoolkit.setting as setting
from qmjob import QMRun

class QMInp(object):
  def __init__(self, 
               structure_inp, 
               program=setting.qmcode, 
               **kwargs):
    self.program = program
    self.atom_count = 0
    if 'info' in kwargs:
      self.info = kwargs['info']
    else: 
      self.info = structure_inp

    # take input 'program' to choose corresponding format
    if self.program == 'cpmd':
      self.inp = cpmd.inp(structure_inp, self.info)
    elif self.program=='vasp':
      self.inp = vasp.inp(structure_inp, self.info)
    else:
      ut.exit('program', self.program, 'is not implemented')

    if 'mode' in kwargs:
      self.inp.setting.mode = kwargs['mode']
      self.inp.setting.set_mode = True
    if 'temperature' in kwargs:
      self.inp.setting.temperature = kwargs['temperature']
    if 'temperature_tolerance' in kwargs:
      self.inp.setting.tolerance =\
      kwargs['temperature_tolerance']
    if 'md_step' in kwargs:
      self.inp.setting.md_step = kwargs['md_step']
    if 'md_sample_rate' in kwargs:
      self.inp.setting.md_sample_rate = kwargs['md_sample_rate']
    if 'theory' in kwargs:
      self.inp.setting.theory = kwargs['theory']
    if 'vdw' in kwargs: 
      self.inp.setting.vdw = kwargs['vdw']
      self.inp.setting.set_vdw = True
    if 'wf_step' in kwargs:
      self.inp.setting.maxstep = kwargs['wf_step']
      self.inp.setting.set_step = True
    if 'cutoff' in kwargs:
      self.inp.setting.cutoff = kwargs['cutoff']
    if 'periodic' in kwargs and kwargs['periodic']:
      self.inp.setting.isolated = False
      self.inp.setting.set_center = False
    if 'celldm' in kwargs:
      self.inp.setting.celldm = kwargs['celldm']
      self.inp.setting.set_celldm = True
    if 'PPext' in kwargs:
      self.inp.setting.PPext = kwargs['PPext']

  def view(self, name=None):
    self.inp.structure.view(name)

  def setAtom(self, atom_list, atom_string):
    self.atom_dict = {}
    if type(atom_list) == int:
      atom_list = [atom_list]
    for i in range(len(atom_list)):
      if atom_list[i] > 0:
        self.atom_dict[i] = atom_string
      else:
        ut.exit("atom index starts from 1")
        
    
#    _tmp = self.atom_count
#    if type(atom_list) == int:
#      atom_list = [atom_list]
#    for I in atom_list:
#      if I == 0:
#        ut.exit("atom index starts from 1")
#      i = I-1
#      self.inp.structure.type_list[i] = atom_string
#      # negative counter for count added atom string
#      self.inp.structure.Z[i] = self.atom_count
#    self.inp.atom_list[str(self.atom_count)] = atom_string
#    self.atom_count = _tmp - 1

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
    if type(structure)==str:
      self.inp.structure = geometry.Molecule()
      self.inp.structure.read(structure)
    else:
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

  def setMode(self, mode, **kwargs):
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

  def setVDW(self, vdw):
    self.inp.setting.vdw = vdw
    self.inp.setting.set_vdw = True

  def setTemperature(self, temperature):
    self.inp.setting.temperature = temperature
  def setTolerance(self, tolerance):
    self.inp.setting.tolerance = tolerance
  def setMDStep(self, md_step):
    self.inp.setting.md_step = md_step
  def setMDSampleRate(self, rate):
    self.inp.setting.md_sample_rate = rate

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

  def run(self, *run_arg, **run_kw):
    return self.inp.run(*run_arg, **run_kw)
