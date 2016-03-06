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
    if 'root_dir' not in kwargs:
      self.setting['root_dir'] = name
    self.setting['reset'] = True

    def writeInp(name=None, **setting):
      setting['no_update'] = True
      inp, molecule = \
        super(PlanewaveInput, self).write(name, **setting)
  
      if 'ks_states' in setting and 'ks_states':
        setting['mode'] = 'ks_states'
 
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
        inp.write(' LANCZOS PARAMETER N=5\n')
        inp.write('  1000  16  20  1.D-9\n')
        inp.write('  0.05          1.D-11\n')
        inp.write('  0.01          1.D-13\n')
        inp.write('  0.0025        1.D-16\n')
        inp.write('  0.001         1.D-18\n')
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
            inp.write(' %9.6f' % component)
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
#      if 'scale' in setting:
#        inp.write(' SCALE SX=%d SY=%d SZ=%d\n' %\
#          (setting['scale'][0],
#           setting['scale'][1],
#           setting['scale'][2]))
      if 'kmesh' in setting:
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
      del setting['save_wf']
      max_wf = max(wf_list)
      if max_wf > self.getValenceElectrons() / 2:
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
    elif 'ks_states' in setting and setting['ks_states']:
      n_ks = setting['ks_states']
      del setting['ks_states']
      setting['no_subfolder'] = False
      setting['save_restart'] = True
      sub_name = name
      if name:
        sub_name = name + '_01'
      inp = writeInp(sub_name, **setting)
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
      root = os.path.split(os.path.abspath(qmout))[0]
      mol_path = os.path.join(root, 'GEOMETRY.xyz')
      try:
        self.molecule = qtk.Molecule(mol_path)
      except:
        pass
      self.getEt(qmout)

  def getEt(self, name):
    out = sys.stdout if re.match('stdout',name)\
          else open(name, 'r')

    done = False
    finished = False
    converged = True
  
    scf_p = re.compile('^ *[0-9]*  [0-9]\.[0-9]{3}E-[0-9]{2}   .*')
    Et_cpmd = re.compile('.*TOTAL ENERGY = *([-0-9\.]*)')
    convergence = re.compile('.*BUT NO CONVERGENCE.*')
    soft_exit = re.compile('.*SOFT EXIT REQUEST.*')
    done_cpmd = re.compile(' \* *FINAL RESULTS *\*')
    qmInfo = re.compile('.*qmInfo.*')
    info_head = re.compile('.*qmInfo:')
    info_tail = re.compile(' *\*\*.*')

    while True:
      line = out.readline()
      if not line: break

      if (re.match(scf_p, line)):
        try:
          data = [float(x) for x in line.split()]
          self.scf_step = int(data[0])
        except:
          qtk.report("\n\nFailed while reading file:", name,
                    'at line: ', line,
                    '... skipping!! \n', color='yellow')
      elif re.match(convergence, line) and self.scf_step > 5:
        converged = False
      elif re.match(soft_exit, line):
        converged = False
      elif re.match(qmInfo, line):
        tmp1 = re.sub(info_head, '', line)
        tmp2 = re.sub(info_tail, '', tmp1)
        self.info = tmp2
      elif (re.match(done_cpmd,line)):
        done = True,
      elif (re.match(Et_cpmd, line)) and done and converged:
        self.Et = float(Et_cpmd.match(line).group(1))
        finished = True

# not used by PP object but by QMInp cpmd parts
def PPString(inp, mol, i, n, outFile):
  """
  append PP file names to inp.pp_files
  """
  alchemy = re.compile('^\w*2\w*_\d\d\d$')
  ppstr = re.sub('\*', '', mol.string[i])
  if ppstr:
    PPStr = '*' + ppstr
    pp_root, pp_ext = os.path.split(ppstr)
    if pp_ext != 'psp':
      PPStr = PPStr + '.psp'
  elif 'vdw' in inp.setting\
  and inp.setting['vdw'].lower == 'dcacp':
    # PPStr: Element_qve_dcacp_theory.psp
    PPStr = '*'+mol.type_list[i] + '_dcacp_' +\
      inp.setting['pp_theory'].lower() + '.psp'
  else:
    # PPStr: Element_qve_theory.psp
    nve = qtk.n2ve(mol.type_list[i])
    PPStr = '*'+mol.type_list[i] + '_q%d_' % nve +\
      inp.setting['pp_theory'].lower() + '.psp'
  outFile.write(PPStr + '\n')
  pp_file_str = re.sub('\*', '', PPStr)
  xc = inp.setting['pp_theory'].lower()
  if not mol.string[i]:
    PPCheck(xc, mol.type_list[i].title(), pp_file_str)
  elif alchemy.match(mol.string[i]):
    alchemyPP(xc, pp_file_str)
  pp_file = os.path.join(qtk.setting.cpmd_pp, pp_file_str)

  if inp.setting['pp_type'].title() == 'Goedecker':
    lmax = 'F'
  outFile.write(' LMAX=%s\n %3d\n' % (lmax, n))
  inp.pp_files.append(re.sub('\*', '', PPStr))

# not used by PP object but by QMInp cpmd parts
def PPCheck(xc, element, pp_file_str, **kwargs):
  if xc == 'lda':
    xc = 'pade'
  ne = qtk.n2ve(element)
  try:
    pp_path = os.path.join(xc, element + '-q' + str(qtk.n2ve(element)))
    pp_file = os.path.join(qtk.setting.cpmd_pp_url, pp_path)
    saved_pp_path = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
    if not os.path.exists(saved_pp_path):
      if pp_file:
        new_pp = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
        pp_content = urllib2.urlopen(pp_file).read()
        qtk.report('', 'pp file %s not found in %s. ' \
                   % (pp_file_str, qtk.setting.cpmd_pp) + \
                   'but found in cp2k page, download now...')
        new_pp_file = open(new_pp, 'w')
        new_pp_file.write(pp_content)
        new_pp_file.close()
        pp_file = new_pp
    return saved_pp_path
  except:
    qtk.warning('something wrong with pseudopotential')

def alchemyPP(xc, pp_file_str):
  pp_path = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
  if not os.path.exists(pp_path):
    root, _ = os.path.splitext(pp_file_str)
    element_str = re.sub('_.*', '', root)
    fraction = float(re.sub('.*_', '', root))/100
    element1 = re.sub('2.*', '', element_str)
    element2 = re.sub('.*2', '', element_str)
    str1 = element1 + "_q" + str(qtk.n2ve(element1)) +\
           "_" + xc + '.psp'
    str2 = element2 + "_q" + str(qtk.n2ve(element2)) +\
           "_" + xc + '.psp'
    pp1 = PP(PPCheck(xc, element1, str1))
    pp2 = PP(PPCheck(xc, element2, str2))
    pp = mutatePP(pp1, pp2, fraction)
    pp.write(pp_path)

def mutatePP(pp1, pp2, fraction):
  if type(pp1) is str:
    if pp1.upper() == 'VOID':
      pp1 = PP()
    else:
      pp1 = PP(pp1)
  if type(pp2) is str:
    if pp2.upper() == 'VOID':
      pp2 = PP()
    else:
      pp2 = PP(pp2)
  pp1 = pp1*(1-fraction)
  pp2 = pp2*fraction
  pp = pp1 + pp2
  if pp1.param['Z']*pp2.param['Z'] > 0:
    if fraction > 0.5:
      pp.param['Z'] = pp2.param['Z']
  else:
    if pp1.param['Z'] == 0:
      pp.param['Z'] = pp2.param['Z']
    else:
      pp.param['Z'] = pp1.param['Z']
  return pp
