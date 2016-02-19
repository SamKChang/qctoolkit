import qctoolkit as qtk
import re, os, shutil, copy, sys
import numpy as np
import qctoolkit.QM.qmjob as qmjob

class InpContent(object):
  """
  general input file content object. It manages checking existance, 
  creating folder, copy necessary files, 
  append extension/prefix/suffix, writing file or print to screen

  Attribute: (default at qtk.setting.file_setup)
    prefix(str),
    suffix(str),
    extension(str),
    root_dir(str),
    path(str),
    output(boolean),
    dependent_files(list str)
    finalized(boolean),
    final_path(str)
    final_name(str)
  """

  def __init__(self, name, **kwargs):
    self.content = []
    self.file_name = name
    self.finalized = False

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
    """
    finalize input content by either 1) creating folder, 
    coping dependent files and write input content to file 
    with correct name or 2) print to screen

    Note that is returns full path to the file is writes to 
    for posprocessing
    """

    if self.output:
      name = self.file_name
      if self.prefix:
        name = self.prefix + name
      if self.suffix:
        name = name + self.suffix 
      if self.extension:
        name = name + '.' + self.extension
      self.final_name = name
      full_dir_path = self.path
      if self.root_dir:
        full_dir_path = os.path.join(self.path, self.root_dir)
      if name:
        full_path = os.path.join(full_dir_path, name)
        if os.path.exists(full_path):
          qtk.prompt(full_path + ' exists, overwrite?')
          try:
            os.remove(full_path)
          except OSError:
            shutil.rmtree(full_path)
          except:
            qtk.exit("can not remove file: " + name)
        if not os.path.exists(full_dir_path):
          os.makedirs(full_dir_path)
        self.final_path = full_dir_path
        for dep in self.dependent_files:
          if os.path.exists(dep):
            shutil.copy(dep, full_dir_path)
          else:
            qtk.warning("dependent file: %s not found" % dep)

    inp = sys.stdout if not self.output else open(full_path, 'w')
    for string in self.content:
      if string:
        inp.write(string)
    if self.output:
      inp.close()
      self.finalized = True

class QMWorker(object):
  def __init__(self, qmcode, **kwargs):
    self.qmcode = qmcode
    self.setting = kwargs

  def start(self, InpClass, new_name=None):
    self.inp = InpClass
    if not new_name:
      self.name = self.inp.final_name
    else:
      self.name = new_name
    if self.inp.finalized:
      cwd = os.getcwd()
      os.chdir(self.inp.final_path)
      run_setup = copy.deepcopy(self.setting)
      del run_setup['program']

      out = qmjob.QMRun(self.name, self.qmcode, **run_setup)
      os.chdir(cwd)
 
#      try:
#        out = qmjob.QMRun(self.name, self.qmcode, **run_setup)
#      except:
#        qtk.warning("qmjob finished unexpectedly for '" + \
#                    self.name + "'")
#        out = GenericQMOutput()
#      finally:
#        os.chdir(cwd)

      return out
    else:
      qtk.exit("InpContent not finalized, no inp name?")

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

  def run(self, qmcode, name=None, **kwargs):
    if not name:
      name = self.molecule.name
    self.setting.update(**kwargs)
    return QMWorker(self.setting['program'], **self.setting), name

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
