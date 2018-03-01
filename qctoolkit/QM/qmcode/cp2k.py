import qctoolkit as qtk
from qctoolkit.QM.planewave_io import PlanewaveInput
from qctoolkit.QM.planewave_io import PlanewaveOutput
import sys, os, re, copy, shutil
import qctoolkit.QM.qmjob as qmjob
from qctoolkit.QM.pseudo.pseudo import PP
import pkg_resources
import numpy as np
import urllib2
import universal as univ
import sys
from urlparse import urlparse
from BeautifulSoup import BeautifulSoup
#reload(sys)
#sys.setdefaultencoding('utf8')

class inp(PlanewaveInput):
  """
  cpmd input class.
  Note!! automatic PP download for DCACPs is not yet done
  but if the PPs are in the right location with right name
  automatic interpolation would work in principle
  """
  __doc__ = PlanewaveInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    PlanewaveInput.__init__(self, molecule, **kwargs)
    self.setting.update(**kwargs)
    if 'ks_states' in kwargs:
      self.setting['mode'] = kwargs['ks_states']
    if 'pp_type' not in kwargs:
      self.setting['pp_type'] = 'Goedecker'
    if 'pp_theory' not in kwargs:
      self.setting['pp_theory'] = self.setting['theory']
    self.backup()

  def run(self, name=None, **kwargs):
    kwargs['no_subfolder'] = False
    if not name:
      kwargs['new_name'] = self.molecule.name
    else:
      kwargs['new_name'] = name
    self.setting.update(kwargs)
    return univ.runCode(self, PlanewaveInput, name, **self.setting)

  def write(self, name=None, **kwargs):
    self.setting.update(kwargs)
    if 'geopt' in self.setting and self.setting['geopt']:
      self.setting['mode'] = 'geopt'
    if 'root_dir' not in kwargs:
      self.setting['root_dir'] = name
    self.setting['reset'] = True

    def writeInp(name=None, **setting):
      setting['no_update'] = True
      inp, molecule = PlanewaveInput.write(self, name, **setting)

      if 'ks_states' in setting and 'ks_states':
        setting['mode'] = 'ks_states'

      if 'big_memory' not in setting:
        if qtk.setting.memory > 16:
          setting['big_memory'] = True

#      if molecule.scale:
#        molecule.R = molecule.R_scale
#        setting['scale'] = molecule.scale
  
      inp.write('&INFO\n')
      inp.write(' ' + setting['info'] + '\n')
      inp.write('&END\n\n')
  
      inp.write('&CPMD\n MIRROR\n')
      if 'save_restart' not in setting\
      or not setting['save_restart']:
        inp.write(' BENCHMARK\n')
        inp.write('  1 0 0 0 0 0 0 0 0 0\n')
      if 'init_random' in setting\
      and setting['init_random']:
        inp.write(' INITIALIZE WAVEFUNCTION RANDOM\n')
      if setting['mode'].lower() == 'single_point':
        inp.write(' OPTIMIZE WAVEFUNCTION\n')
      elif setting['mode'].lower() == 'geopt':
        inp.write(' OPTIMIZE GEOMETRY\n')
      elif setting['mode'].lower() == 'ks_states':
        if not setting['restart']:
          setting['restart'] = True
        inp.write(' KOHN-SHAM ENERGIES\n  %d\n'\
          % setting['ks_states'])
        #inp.write(' LANCZOS PARAMETER N=5\n')
        #inp.write('  1000  16  20  1.D-9\n')
        #inp.write('  0.05          1.D-11\n')
        #inp.write('  0.01          1.D-13\n')
        #inp.write('  0.0025        1.D-16\n')
        #inp.write('  0.001         1.D-18\n')
      elif setting['mode'].lower() == 'md':
        inp.write(' MOLECULAR DYNAMICS BO\n')
      if re.match(setting['mode'], 'md'):
        inp.write(' TRAJECTORY SAMPLE XYZ\n  %d\n'\
          % setting['sample_period'])
        inp.write(' TEMPERATURE\n')
        inp.write('  %.1f\n' % setting['T'])
        if setting['thermostat'].lower() == 'nose-hoover':
          inp.write(' NOSE IONS\n  %.1f %.1f\n' % \
            (setting['T'], setting['T_tolerance']))
      if 'geometry_convergence' in setting:
        inp.write(' CONVERGENCE GEOMETRY\n  %.2E\n' %\
          setting['geometry_convergence'])
      if 'md_step' in setting:
        inp.write(' MAXSTEP\n  %d\n' % setting['md_step'])
      if 'wf_convergence' in setting:
        inp.write(' CONVERGENCE ORBITAL\n  %.2E\n' %\
          setting['wf_convergence'])
      if 'scf_step' in setting:
        inp.write(' MAXITER\n  %d\n' % setting['scf_step'])
      if 'restart' in setting and setting['restart']:
        inp.write(' RESTART WAVEFUNCTION\n')
      if 'restart_wavefunction_file' in setting:
        if os.path.exists(setting['restart_wavefunction_file']):
          inp.dependent_files.append(
            [setting['restart_wavefunction_file'], 'RESTART']
          )
        else:
          qtk.warning('%s not found' % setting['restart_wavefunction_file'])
      if 'fix_molecule' in setting\
       and setting['fix_molecule']:
        inp.write(' CENTER MOLECULE OFF\n')
      if ('save_density' in setting and setting['save_density'])\
      or ('save_wf' in setting and setting['save_wf']):
        inp.write(' RHOOUT')
        if setting['save_wf']:
          inp.write(' BANDS\n  %d\n ' % len(setting['save_wf']))
          for i in setting['save_wf']:
            inp.write(' %d' % i)
        inp.write('\n')

      if molecule.multiplicity != 1:
        inp.write(' LOCAL SPIN DENSITY\n') 
      if 'big_memory' in setting and setting['big_memory']:
        inp.write(' MEMORY BIG\n')
      inp.write('&END\n\n')
  
      inp.write('&DFT\n FUNCTIONAL %s\n&END\n\n' % \
        setting['theory'].upper())
  
      inp.write('&SYSTEM\n')
      if setting['periodic']:
        inp.write(' CELL VECTORS\n')
        lattice_vec = self.setting['lattice']
        for vec in lattice_vec:
          inp.write(' ')
          for component in vec:
            inp.write(' % 11.6f' % component)
          inp.write('\n')
      else:
        inp.write(' SYMMETRY\n')
        inp.write('  ISOLATED\n')
        if setting['isolation'] == 'mt':
          inp.write(' POISSON SOLVER TUCKERMAN\n')
        inp.write(' CELL ABSOLUTE\n ')
        for d in setting['celldm']:
          inp.write(' %6.3f'% float(d))
        inp.write('\n')
      if setting['unit'].lower() == 'angstrom':
        inp.write(' ANGSTROM\n')
      inp.write(' CUTOFF\n  %.1f\n' % setting['cutoff'])
      if 'kmesh' in setting and setting['kmesh']:
        inp.write(' KPOINTS MONKHORST-PACK\n  %d %d %d\n' %\
          (setting['kmesh'][0],
           setting['kmesh'][1],
           setting['kmesh'][2]))
      if 'mesh' in setting:
        inp.write(' MESH\n  %d %d %d\n' %\
          (setting['mesh'][0],
           setting['mesh'][1],
           setting['mesh'][2]))
      if molecule.charge != 0:
        inp.write(' CHARGE\n  %d\n' % molecule.charge)
      if molecule.multiplicity != 1:
        inp.write(' MULTIPLICITY\n  %d\n' % molecule.multiplicity)
      inp.write('&END\n\n')
  
      inp.write('&ATOMS\n')
      molecule.sort()
      type_index = molecule.index
      type_list = molecule.type_list
      # loop through sorted atom type indices
      for atom_type in xrange(0,len(type_index)-1):
        # number of each atom type
        type_n = type_index[atom_type+1] - type_index[atom_type]
        PPString(self, molecule, type_index[atom_type], type_n, inp)
        for I in\
        xrange(type_index[atom_type],type_index[atom_type+1]):
          inp.write('  ')
          for i in range(3):
            inp.write(' % 8.4f' % molecule.R[I,i])
          inp.write('\n')
      inp.write('&END\n')

      for pp in self.pp_files:
        pp_file = os.path.join(qtk.setting.cpmd_pp, pp)
        inp.dependent_files.append(pp_file)
  
      if 'no_cleanup' in setting and setting['no_cleanup']:
        inp.close(no_cleanup=True)
      else:
        inp.close()

      return inp
  
      #return univ.writeReturn(inp, name, **setting)

    inp = None
    setting = copy.deepcopy(self.setting)
    if 'save_wf' in setting and setting['save_wf']:
      if type(setting['save_wf']) is int:
        setting['save_wf'] = [setting['save_wf']]
      wf_list = setting['save_wf']
      max_wf = max(wf_list)
      if max_wf > self.getValenceElectrons() / 2:
        del setting['save_wf']
        setting['no_subfolder'] = False
        setting['save_restart'] = True
        sub_name = name
        if name:
          sub_name = name + '_01'
        inp = writeInp(sub_name, **setting)
        setting['ks_states'] = max_wf - self.getValenceElectrons() / 2
        setting['no_cleanup'] = True
        setting['save_wf'] = wf_list
        if name:
          sub_name = name + '_02'
        writeInp(sub_name, **setting)
      else:
        inp = writeInp(name, **setting)
    elif 'ks_states' in setting and setting['ks_states']\
    and not setting['restart']:
      n_ks = setting['ks_states']
      del setting['ks_states']
      setting['mode'] = 'single_point'
      setting['no_subfolder'] = False
      setting['save_restart'] = True
      sub_name = name
      if name:
        sub_name = name + '_01'
      inp = writeInp(sub_name, **setting)
      setting['mode'] = 'ks_states'
      setting['ks_states'] = n_ks
      setting['no_cleanup'] = True
      if name:
        sub_name = name + '_02'
      writeInp(sub_name, **setting)
    else:
      inp = writeInp(name, **setting)

    return inp

class out(PlanewaveOutput):
  def __init__(self, qmout, **kwargs):
    PlanewaveOutput.__init__(self, qmout, **kwargs)
    self.info = ''
    if qmout:
      out = open(qmout, 'r')
      data = out.readlines()
      out.close()

      self.details = filter(lambda x: 'energy' in x and ':' in x, data)

      Et = filter(lambda x: 'Total energy:' in x, self.details)[0]
      self.Et = float(Et.split()[-1])
      self.energies = {}
      for term in self.details:
        words = term.split()
        k, v = ' '.join(words[:-1]).replace(':', ''), float(words[-1])
        self.energies[k] = v
