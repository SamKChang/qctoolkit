import qctoolkit as qtk
import re, os, shutil, copy
import numpy as np

class GenericQMInput(object):
  def __init__(self, molecule, **kwargs):
    self.setting = kwargs
    if 'program' not in kwargs:
      self.setting['program'] = qtk.setting.qmcode
    else:
      self.setting['program'] = kwargs['program']
    self.molecule = qtk.Structure(molecule)

    if 'info' not in kwargs:
      self.setting['info'] = self.molecule.name
    if 'theory' not in kwargs:
      self.setting['theory'] = 'pbe'
    if 'mode' not in kwargs:
      self.setting['mode'] = 'single_point'
    if 'scf_step' not in kwargs:
      self.setting['scf_step'] = 1000
    if 'wf_convergence' not in kwargs:
      self.setting['wf_convergence'] = 1E-5
    if 'fix_molecule' in kwargs:
      self.setting['fix_molecule'] = kwargs['fix_molecule']
    if 'unit' not in kwargs:
      self.setting['unit'] = 'angstrom'
    if re.match(self.setting['mode'].lower(), 'md'):
      if 'T' not in kwargs:
        self.setting['T'] = 298
      if 'thermostat' not in kwargs:
        self.setting['thermostat'] = 'Nose-Hoover'
      if 'T_tolerance' not in kwargs:
        self.setting['T_tolerance'] = 0.1
      if 'sample_period' not in kwargs:
        self.setting['sample_period'] = 10
      if 'md_step' not in kwargs:
        self.setting['md_step'] = 1000

    if self.setting['mode'] == 'geopt':
      self.setting['geometry_convergence'] = 1E-4

    if 'save_density' in kwargs:
      self.setting['save_density'] = kwargs['save_density']
    if 'save_restart' in kwargs:
      self.setting['save_restart'] = kwargs['save_restart']
    if 'restart' in kwargs:
      self.setting['restart'] = kwargs['restart']

    if 'prefix' in kwargs:
      self.setting['prefix'] = kwargs['prefix']
    if 'suffix' in kwargs:
      self.setting['suffix'] = kwargs['suffix']

  def __repr__(self):
    return self.molecule.name + ': ' + self.setting['program']

  def view(self, name=None):
    self.molecule.view(name)

  def setAtom(self, *args, **kwargs):
    self.molecule.setAtom(*args, **kwargs)

  def removeAtom(self, index):
    self.molecule.remove_atom(index)

  def isolateAtoms(self, indices):
    self.molecule.isolate_atoms(indices)

  def setChargeMultiplicity(self, *args, **kwargs):
    self.molecule.setChargeMultiplicity(*args, **kwargs)

  def create_folder(self, name=None):
    if not name:
      name = self.molecule.name
    if os.path.exists(name) and not qtk.setting.no_warning:
      qtk.prompt(name + ' exists, overwrite?')
      shutil.rmtree(name)
    os.mkdir(name)
    return name

  def cm_check(self, mol):
    ve = mol.getValenceElectrons()
    if (ve % 2) == (mol.multiplicity % 2):
      msg = "Multiplicity %d " % mol.multiplicity + \
            "and %d valence electrons " % ve +\
            "\n(with charge %3.1f) " % float(mol.charge) +\
            "are not compatible"
      qtk.exit(msg)

  def run(self):
    raise NotImplementedError("Please Implement run method")

  def write(self):
    raise NotImplementedError("Please Implement write method")

class GenericQMOutput(object):
  def __init__(self, output, **kwargs):
    self.Et = np.nan
    self.scf_step = np.nan
    self.unit = 'Eh'

  def __repr__(self):
    return str(self.Et)

  def inUnit(self, unit):
    unitStr = self.unit + '-' + unit
    if not unitStr.lower() == 'eh-eh':
      new_E = qtk.convE(self.Et, unitStr)
      self.Et = new_E
      return new_E
    else: return self.Et

  def __add__(self, other):
    out = copy.deepcopy(self)
    if isinstance(other, qtk.QM.general_io.GenericQMOutput):
      if self.unit == other.unit:
        out.Et = self.Et + other.Et
      else:
        unitStr = self.unit + '-' + other.unit
        out.Et = self.Et + qtk.convE(other.Et, unitStr)
      out.scp_step = max(self.scf_step, other.scf_step)
    elif type(other) is int or float:
      out.Et = self.Et + other
    return out

  def __sub__(self, other):
    out = copy.deepcopy(self)
    if isinstance(other, qtk.QM.general_io.GenericQMOutput):
      if self.unit == other.unit:
        out.Et = self.Et - other.Et
      else:
        unitStr = self.unit + '-' + other.unit
        out.Et = self.Et + qtk.convE(other.Et, unitStr)
      out.scp_step = max(self.scf_step, other.scf_step)
    elif type(other) is int or float:
      out.Et = self.Et - other
    return out

  def __mul__(self, other):
    out = copy.deepcopy(self)
    if type(other) is int or float:
      out.Et = self.Et * other
    return out

  def __div__(self, other):
    out = copy.deepcopy(self)
    if type(other) is int or float:
      out.Et = self.Et / other
    return out
