import qctoolkit as qtk
import re, os, shutil, copy, sys
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import universal as univ

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
    self.reset()
  """

  def __init__(self, name, **kwargs):
    self.content = []
    self.file_name = name

    setup = copy.deepcopy(qtk.setting.file_setup)

    for string, value in setup.iteritems():
      if string in kwargs:
        setattr(self, string, kwargs[string])
      else:
        setattr(self, string, value)

    if 'dependent_files' in kwargs:
      self.dependent_files = kwargs['dependent_files']
      if type(self.dependent_files) is not list:
        self.dependent_files = [self.dependent_files]
    if 'output' not in kwargs:
      if self.file_name:
        self.output = True
    if 'root_dir' in kwargs:
      self.root_dir = kwargs['root_dir']
    self.path = os.getcwd()

  def write(self, string):
    self.content.append(string)

  def cleanPath(self, path):
    if os.path.exists(path):
      qtk.prompt(path + ' exists, overwrite?')
      try:
        os.remove(path)
      except OSError:
        shutil.rmtree(path)
      except:
        qtk.warning("can not remove file: " + path)
    
  def close(self, **kwargs):
    """
    finalize input content by either 1) creating folder, 
    coping dependent files and write input content to file 
    with correct name or 2) print to screen

    Note that is returns full path to the file is writes to 
    for posprocessing
    """

    if self.output:

      # process file name
      name = self.file_name
      if self.prefix:
        name = self.prefix + name
      if self.suffix:
        name = name + self.suffix 
      if self.extension:
        name = name + '.' + self.extension
      self.final_name = name
      full_dir_path = self.path

      # process directory name and clean up
      if self.root_dir:
        full_dir_path = os.path.join(self.path, self.root_dir)
      if name:
        full_path = os.path.join(full_dir_path, name)
        if not 'no_cleanup' in kwargs or not kwargs['no_cleanup']:
          if self.root_dir:
            self.cleanPath(full_dir_path)
        self.cleanPath(full_path)
        if not os.path.exists(full_dir_path):
          os.makedirs(full_dir_path)
        self.final_path = full_dir_path

        # copy dependent files
        if 'dependent_files' in kwargs:
          self.dependent_files.extend(kwargs['dependent_files'])
        for dep in set(self.dependent_files):
          if os.path.exists(dep):
            shutil.copy(dep, full_dir_path)
          else:
            qtk.warning('dependent file: %s not found' % dep)

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

      if 'debug' in self.setting and self.setting['debug']:
        out = qmjob.QMRun(self.name, self.qmcode, **run_setup)
        os.chdir(cwd)
      else:
        try:
          out = qmjob.QMRun(self.name, self.qmcode, **run_setup)
        except:
          qtk.warning("qmjob finished unexpectedly for '" + \
                      self.name + "'")
          out = GenericQMOutput()
        finally:
          os.chdir(cwd)

      return out
    else:
      qtk.exit("InpContent not finalized, no inp name?")

class GenericQMInput(object):
  """
  From GenericQMInput:
  holder class of general QM job. It initiates by Molecule class and
  the choice of qmcode (default at qtk.setting.qmcode). 

  attributes:
    setting(dict) --- contain all qm setup
    molecule(Molecule) --- for geometry, atom types, charges

  methods:
    write(name) --- write input file into a file/foler called 'name'
                    extension, prefix, suffix can be set by kwargs
    run(name) --- create a folder call 'name', go into that folder,
                  write an input file, start qmjob
  ===
  """
  def __init__(self, molecule, **kwargs):

    self.setting = copy.deepcopy(kwargs)
    self.molecule = copy.deepcopy(molecule)

    # local variable to construct 'setting' for parsing
    setup = copy.deepcopy(qtk.setting.qm_setup)
    md_setup = copy.deepcopy(qtk.setting.md_setup)

    for string, value in setup.iteritems():
      if string in kwargs:
        self.setting[string] = kwargs[string]
      else:
        self.setting[string] = value
    if self.molecule.periodic:
      if 'periodic' not in kwargs:
        self.setting['periodic'] = True

    if self.setting['mode'] == 'md':
      for string, value in md_setup.iteritems():
        if string in kwargs:
          self.setting[string] = kwargs[string]
        else:
          self.setting[string] = value

    if self.setting['mode'] != 'geopt':
      del self.setting['geometry_convergence']
    self.setting['info'] = self.molecule.name

    method_list = [
                    method for method in dir(molecule)\
                    if callable(getattr(molecule, method))
                  ]
    method_list = filter(lambda x: '__' not in x, method_list)
    method_list = filter(lambda x: 'write' not in x, method_list)
    for method in method_list:
      setattr(self, method, getattr(self.molecule, method))

  def reset(self):
    if self.setting_backup:
      self.setting = copy.deepcopy(self.setting_backup)
      self.molecule = copy.deepcopy(self.molecule_backup)

  def backup(self):
    self.setting_backup = copy.deepcopy(self.setting)
    self.molecule_backup = copy.deepcopy(self.molecule)

  def __repr__(self):
    return self.molecule.name + ': ' + self.setting['program']

  def cm_check(self, mol):
    ve = mol.getValenceElectrons()
    if (ve % 2) == (mol.multiplicity % 2):
      msg = "Multiplicity %d " % mol.multiplicity + \
            "and %d valence electrons " % ve +\
            "\n(with charge %3.1f) " % float(mol.charge) +\
            "are not compatible"
      qtk.exit(msg)

  def run(self, qmcode, name=None, **kwargs):
    self.setting.update(kwargs)
    if not name:
      name = self.molecule.name
    return QMWorker(self.setting['program'], **self.setting), name

  def write(self, name, **kwargs):

    if 'no_reset' in kwargs and kwargs['no_reset']:
      self.reset()
    kwargs.update(self.setting)

    # unify output name/directory
    if name:
      name = os.path.splitext(name)[0]

    # unify periodicity
    cell_check = False
    if 'celldm' in self.setting:
      cell_check = copy.deepcopy(self.setting['celldm'])
      if type(cell_check) is np.ndarray:
        cell_check = True
    if not self.setting['periodic']:
      prop_list = ['celldm', 'scale', 'box']
      for prop in prop_list:
        if prop in self.setting:
          self.setting[prop] = False
    elif not cell_check:
      self.setting['celldm'] = self.setCelldm()
      self.setting['box'] = self.setting['celldm'][:3]

    # unify wavefunction output
    if 'save_wf' in self.setting\
    and type(self.setting['save_wf']) is int:
      self.setting['save_wf'] = [self.setting['save_wf']]

    # set default root_dir ONLY if it is not set
    setting = copy.deepcopy(self.setting)
    if 'root_dir' not in kwargs:
      if name:
        setting['root_dir'] = name
      else:
        setting['root_dir'] = self.molecule.name
    else:
      setting['root_dir'] = kwargs['root_dir']
      del kwargs['root_dir']

    if 'no_update' not in kwargs or not kwargs['no_update']:
      self.setting.update(kwargs)

    inp = InpContent(name, **setting)
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
      self.Et, self.unit = qtk.convE(self.Et, unitStr, '-')
    return self.Et

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
