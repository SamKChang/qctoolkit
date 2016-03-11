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
  espresso input class.
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

    mode_dict = {
                  'single_point': 'scf'
                }

    self.molecule.sort()
    type_index = self.molecule.index
    type_list = self.molecule.type_list

    def writeInp(name=None, **setting):
      inp, molecule = \
        super(PlanewaveInput, self).write(name, **setting)
  
      mode = mode_dict[setting['mode']]

      inp.write("&control\n")
      inp.write(" calculation = '%s',\n" % mode)
      inp.write(" pseudo_dir = './',\n")
      if setting['restart']:
        pass
      else:
        inp.write(" restart_mode = 'from_scratch',\n")
      inp.write(" outdir = 'out',\n")
      inp.write("/\n")

      inp.write("&system\n")
      if not setting['periodic']:
        inp.write(" assume_isolated = 'mt',\n")
      inp.write(" ibrav = 0,\n")
      inp.write(" nat = %d,\n" % molecule.N)
      inp.write(" ntyp = %d,\n" % (len(type_index) - 1))
      inp.write(" ecutwfc = %.1f,\n" % setting['cutoff'])
      inp.write(" nosym = .true.,\n")
      inp.write(" noinv = .true.,\n")
      if molecule.charge != 0:
        inp.write(" tot_charge = %.2f,\n" % molecule.charge)
      inp.write("/\n")

      inp.write("&electrons\n")
      inp.write(" mixing_beta = 0.7,\n")
      inp.write(" diagonalization = 'david',\n")
      inp.write(" conv_thr =  1.0d-9,\n")
      inp.write("/\n")

      inp.write("ATOMIC_SPECIES\n")
      for a in range(len(type_index)-1):
        inp.write(' %2s %8.3f %s\n' % \
          (type_list[type_index[a]], 1.0, 'pp_file'))
      inp.write("\n")

      inp.write("ATOMIC_POSITIONS angstrom\n")
      for a in range(len(type_list)):
        inp.write(' %2s' % type_list[a])
        for i in range(3):
          inp.write(' % 12.8f' % molecule.R[a, i])
        inp.write("\n")
      inp.write("\n")

      inp.write("K_POINTS gamma\n")
      inp.write("\n")

      inp.write("CELL_PARAMETERS angstrom\n")
      lattice_vec = self.setting['lattice']
      for vec in lattice_vec:
        for component in vec:
          inp.write(' %9.6f' % component)
        inp.write('\n')
  
      for pp in self.pp_files:
        pp_file = os.path.join(qtk.setting.cpmd_pp, pp)
        inp.dependent_files.append(pp_file)

      inp.close()

    setting = copy.deepcopy(self.setting)
    writeInp(name, **setting)

class out(PlanewaveOutput):
  def __init__(self, qmout, **kwargs):
    pass

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
    element = mol.type_list[i]
    if 'valence_electrons' in inp.setting\
    and element in inp.setting['valence_electrons']:
      nve = inp.setting['valence_electrons'][element]
    else:
      nve = qtk.n2ve(element)
      
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
