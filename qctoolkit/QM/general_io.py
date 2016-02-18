import qctoolkit as qtk
import re, os, shutil, copy, sys
import numpy as np

class InpContent(object):
  def __init__(self, name, **kwargs):
    self.content = []
    self.file_name = name

    setup = copy.deepcopy(qtk.setting.file_setup)

    for string, value in setup.iteritems():
      if string in kwargs:
        setattr(self, string, kwargs[string])
      else:
        setattr(self, string, value)

    if 'output' not in kwargs:
      if self.file_name or self.root_dir:
        self.output = True

  def write(self, string):
    self.content.append(string)

  def close(self):
    if self.output:
      name = self.file_name
      if self.prefix:
        name = self.prefix + name
      if self.suffix:
        name = name + self.suffix 
      if self.extension:
        name = name + '.' + self.extension
      full_dir_path = self.path
      if self.root_dir:
        full_dir_path = os.path.join(self.path, self.root_dir)
      if name:
        full_path = os.path.join(full_dir_path, name)
        if not os.path.exists(full_dir_path):
          os.makedirs(full_dir_path)
        if os.path.exists(full_path):
          qtk.prompt(name + ' exists, overwrite?')
          try:
            os.remove(name)
          except:
            qtk.exit("can not remove file: " + name)

    inp = sys.stdout if not self.output else open(full_path, 'w')
    for string in self.content:
      if string:
        inp.write(string)
    if self.output:
      inp.close()

class GenericQMInput(object):
  def __init__(self, molecule, **kwargs):
    self.setting = kwargs
    self.molecule = molecule

    # local variable to construct 'setting' for parsing
    setup = copy.deepcopy(qtk.setting.qm_setup)
    md_setup = copy.deepcopy(qtk.setting.md_setup)

    for string, value in setup.iteritems():
      if string in kwargs:
        self.setting[string] = kwargs[string]
      else:
        self.setting[string] = value

    if self.setting['mode'] == 'md':
      for string, value in md_setup.iteritems():
        if string in kwargs:
          self.setting[string] = kwargs[string]
        else:
          self.setting[string] = value

    if self.setting['mode'] != 'geopt':
      del self.setting['geometry_convergence']
    self.setting['info'] = self.molecule.name

  def __repr__(self):
    return self.molecule.name + ': ' + self.setting['program']

  # interface to Molecule class
  def view(self, name=None):
    self.molecule.view(name)

  def setAtoms(self, *args, **kwargs):
    self.molecule.setAtoms(*args, **kwargs)

  def removeAtoms(self, index):
    self.molecule.removeAtoms(index)

  def isolateAtoms(self, indices):
    self.molecule.isolateAtoms(indices)

  def setChargeMultiplicity(self, *args, **kwargs):
    self.molecule.setChargeMultiplicity(*args, **kwargs)
  # end of interface

  def create_folder(self, name=None):
    if not name:
      name = self.molecule.name
    if os.path.exists(name) and not qtk.setting.no_warning:
      qtk.prompt(name + ' exists, overwrite?')
      shutil.rmtree(name)
    os.mkdir(name)
    if 'restart_file' in self.setting:
      shutil.copy(self.setting['restart_file'], name)
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

  def write(self, name, **kwargs):
    self.setting.update(kwargs)
    if name:
      name = os.path.splitext(name)[0]
    inp = InpContent(name, **self.setting)
    molecule = copy.deepcopy(self.molecule)
    self.cm_check(molecule)
    if 'no_molecule' in kwargs and kwargs['no_molecule']:
      return inp
    else:
      return inp, molecule

class GenericQMOutput(object):
  def __init__(self, output=None, **kwargs):
    self.Et = np.nan
    self.scf_step = np.nan
    self.unit = 'Eh'
    if output:
      self.path = qtk.getPath(output)
      _file = qtk.fileStrip(output)
      self.stem, self.ext = os.path.splitext(_file)

  def __repr__(self):
    return str(self.Et)

  def inUnit(self, unit):
    unitStr = self.unit + '-' + unit
    if not unitStr.lower() == 'eh-eh':
      new_E = qtk.convE(self.Et, unitStr, '-')
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
        out.Et = self.Et + qtk.convE(other.Et, unitStr, '-')
      out.scp_step = max(self.scf_step, other.scf_step)
    elif type(other) is (int or float):
      out.Et = self.Et + other
    else:
      out.Et = np.nan
    return out

  def __sub__(self, other):
    out = copy.deepcopy(self)
    if isinstance(other, qtk.QM.general_io.GenericQMOutput):
      if self.unit == other.unit:
        out.Et = self.Et - other.Et
      else:
        unitStr = self.unit + '-' + other.unit
        out.Et = self.Et + qtk.convE(other.Et, unitStr, '-')
      out.scp_step = max(self.scf_step, other.scf_step)
    elif type(other) is (int or float):
      print type(other)
      out.Et = self.Et - other
    else:
      out.Et = np.nan
    return out

  def __mul__(self, other):
    out = copy.deepcopy(self)
    if type(other) is (int or float):
      out.Et = self.Et * other
    else:
      out.Et = np.nan
    return out

  def __div__(self, other):
    out = copy.deepcopy(self)
    if type(other) is (int or float):
      out.Et = self.Et / other
    else:
      out.Et = np.nan
    return out
