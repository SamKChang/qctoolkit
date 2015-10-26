import qctoolkit as qtk
import re, os, shutil

class GenericInput(object):
  def __init__(self, molecule, **kwargs):
    self.setting = {}
    if 'program' not in kwargs:
      self.setting['program'] = qtk.setting.qmcode
    else:
      self.setting['program'] = kwargs['program']
    self.molecule = qtk.Structure(molecule)

    self.setting = kwargs

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
    ve = mol.getValenceElectrons() - mol.charge
    if (ve % 2) == (mol.multiplicity % 2):
      msg = "Multiplicity %d " % self.multiplicity + \
            "and %d valence electrons " % nve +\
            "\n(with charge %3.1f) " % float(self.charge) +\
            "are not compatible"
      qtk.exit(msg)

  def run(self):
    raise NotImplementedError("Please Implement run method")

  def write(self):
    raise NotImplementedError("Please Implement write method")
