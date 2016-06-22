import qctoolkit as qtk
from qctoolkit.QM.planewave_io import PlanewaveInput
from qctoolkit.QM.planewave_io import PlanewaveOutput
import sys, os, re, copy, shutil
import qctoolkit.QM.qmjob as qmjob
from qctoolkit.QM.pseudo.pseudo import PP
import pkg_resources
from collections import OrderedDict as odict
import numpy as np
import urllib2
import glob
import universal as univ
import xml.etree.ElementTree as ET
from cpmd import alchemyPP
from cpmd import PPName
from cpmd import PPCheck as PPCheck_cpmd
import subprocess as sp
#from cpmd import PPCheck as cpCheck
#from cpmd import mutatePP

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
    if 'pp_theory' not in kwargs:
      self.setting['pp_theory'] = self.setting['theory']
    self.backup()

    mode_dict = {
      'single_point': 'scf',
    }

    self.content = odict()

    mode = mode_dict[self.setting['mode']]
    self.content['control'] = odict([
      ('calculation', mode),
      ('pseudo_dir', './'),
    ])
    self.content['system'] = odict([
      ('ibrav', 0),
      ('ecutwfc', self.setting['cutoff']),
    ])
    self.content['electrons'] = odict([
      ('electron_maxstep', self.setting['scf_step']),
    ])

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

    dft_dict = {
      'hse06': 'hse',
      'lda': 'pz',
    }

    def atomIndex(molecule):
      for i in range(len(molecule.index) - 1):
        start = molecule.index[i]
        end = molecule.index[i+1]
        if end - start > 1:
          for I in range(start, end):
            element = molecule.type_list[I] + str(I)
            molecule.setAtoms(I, element = element)
      molecule.sort()
      return molecule
      

    def refCopy(name=None, **setting):

      inp, molecule = \
        PlanewaveInput.write(self, name, **setting)
      molecule.sort()
      molecule = atomIndex(molecule)
      type_index = molecule.index

      def fileInsert(flag, string, data, before=False):
        matched = filter(lambda x: flag in x, data)[0]
        ind = data.index(matched)
        if before:
          data.insert(ind, string)
        else:
          data.insert(ind+1, string)

      def fileModify(flag, value, data):
        if type(value) is not str:
          value = str(value) + ","
        matched = filter(lambda x: flag in x, data)[0]
        ind = data.index(matched)
        str_list = matched.split(' ')
        str_list[-1] = value
        data[ind] = ' '.join(str_list) + '\n'

      def fileGrep(flag, data):
        matched = filter(lambda x: flag in x, data)[0]
        matched = matched.split('=')[1]
        return matched.split(',')[0]

      ref = setting['ref_dir']
      ref_root = os.path.split(ref)[1]
      files = glob.glob(ref + '/*')
      restart = filter(lambda x: '.restart' in x, files)
      rst_inp = os.path.join(ref, ref_root) + '.inp'
      rst_wfc = filter(lambda x: '.wfc' in x, files)
      rst_save = filter(lambda x: '.save' in x, files)[0]
      tar_save = os.path.split(rst_save)[1]

      if len(restart) == 0:
        qtk.exit('no espresso restart file found')
      for rst in restart:
        inp.dependent_files.append(rst)
      for wfc in rst_wfc:
        inp.dependent_files.append(wfc)
      inp.dependent_files.append(rst_save)

      inp_file = open(rst_inp)
      inp_data = inp_file.readlines()
      inp_file.close()
      nat = int(fileGrep('nat =', inp_data))
      ntyp = int(fileGrep('ntyp =', inp_data))

      tar_pp_files = []
      for a in range(len(type_index)-1):
        type_n = type_index[a+1] - type_index[a]
        PPStr = PPString(self, molecule,
          type_index[a], type_n, inp)
        tar_pp_files.append(PPStr)

      species_str = filter(lambda x: 'ATOMIC_SPECIES' in x, inp_data)[0]
      ind = inp_data.index(species_str)
      ref_pp_files = []
      for a in range(ind+1, ind+ntyp+1):
        PPStr = inp_data[a].split(' ')[-1].replace('\n', '')
        ref_pp_files.append(PPStr)

      # add nat, ntyp check
      if molecule.N != nat or len(type_index)-1 != ntyp:
        qtk.exit("number of atoms and atom types must matched")

      for i in range(len(tar_pp_files)):
        tar_pp = tar_pp_files[i]
        ref_pp = ref_pp_files[i]
        pp_file = os.path.join(qtk.setting.espresso_pp, tar_pp)
        if tar_pp not in inp.dependent_files:
          name = os.path.join(tar_save, ref_pp)
          inp.dependent_files.append({pp_file: name})

      try:
        fileModify('electron_maxstep', 1, inp_data)
      except:
        fileInsert('&electrons', ' electron_maxstep = 1,\n', inp_data)
      try:
        fileModify('restart_mode', "'restart',", inp_data)
      except:
        fileInsert('&control', " restart_mode = 'restart',\n", inp_data)

      try:
        fileModify('ALCHEMY', "prediction", inp_data)
      except:
        inp_data.append('\nALCHEMY prediction\n')

      for s in inp_data:
        inp.write(s)
      inp.close()
      return inp

    def writeInp(name=None, **setting):

      inp, molecule = \
        PlanewaveInput.write(self, name, **setting)

      molecule.sort()
      if 'save_restart' in setting and setting['save_restart']:
        molecule = atomIndex(molecule)

      type_index = molecule.index
      type_list = molecule.type_list
      pp_files = []

      self.content['system']['nat'] = molecule.N
      self.content['system']['ntyp'] = len(type_index) - 1
      if setting['full_kmesh']:
        self.content['system']['nosym'] = True
        self.content['system']['noinv'] = True
  
      if 'restart' in setting and setting['restart']:
        self.content['control']['restart_mode'] = 'restart'

      if 'save_density' in setting and setting['save_density']:
        self.content['control']['wf_collect'] = True

      if 'save_wf' in setting and setting['save_wf']:
        self.content['control']['wf_collect'] = True

      if not setting['periodic']:
        self.content['system']['assume_isolated'] = 'mt'

      if 'exx' in setting and setting['exx'] == 'anisotropic':
        self.content['system']['exxdiv_treatment'] = 'vcut_ws'
        self.content['system']['ecutvcut'] = 0.7
        self.content['system']['x_gamma_extrapolation'] = False

      if molecule.charge != 0:
        self.content['system']['tot_charge'] = molecule.charge

      if 'theory' in setting:
        if self.setting['theory'] in dft_dict:
          self.content['system']['input_dft'] = \
            dft_dict[self.setting['theory']]
        else:
          self.content['system']['input_dft'] = self.setting['theory']

      if 'scf_step' in setting:
        self.content['electrons']['electron_maxstep'] =\
          setting['scf_step']

      if 'ks_states' in setting and setting['ks_states']:
        vs = int(round(self.molecule.getValenceElectrons() / 2.0))
        self.content['system']['nbnd'] = setting['ks_states'] + vs
        if 'd_shell' in setting:
          for a in molecule.type_list:
            if a in setting['d_shell'] and qtk.n2ve(a) < 10:
              self.content['system']['nbnd'] += 5

      if 'symmetry' in setting:
        if setting['symmetry'] == 'fcc':
          self.content['system']['ibrav'] = 2
          dm = ['A', 'B', 'C', 'cosBC', 'cosAC', 'cosAB']
          for i in range(6):
            self.content['system'][dm[i]] = float(self
              .molecule.celldm[i])

      for section_key in self.content.iterkeys():
        section = '&' + section_key + '\n'
        inp.write(section)
        for key, value in self.content[section_key].iteritems():
          if type(value) is str:
            entry = " %s = '%s',\n" % (key, value)
          elif type(value) is int:
            entry = ' %s = %d,\n' % (key, value)
          elif type(value) is float:
            entry = ' %s = %14.8E,\n' % (key, value)
          elif type(value) is bool:
            if value:
              entry = ' %s = .true.,\n' % key
            else:
              entry = ' %s = .false.,\n' % key
          inp.write(entry)
        inp.write('/\n')

      inp.write("ATOMIC_SPECIES\n")
      for a in range(len(type_index)-1):
        type_n = type_index[a+1] - type_index[a]
        PPStr = PPString(self, molecule, 
          type_index[a], type_n, inp)
        stem, ext = os.path.splitext(PPStr)
        if ext != '.UPF':
          PPStr = PPStr + '.UPF'
        mass = qtk.n2m(type_list[type_index[a]])
        inp.write(' %-3s % 6.3f %s\n' % \
          (type_list[type_index[a]], mass, PPStr))
        pp_files.append(PPStr)
      inp.write("\n")

      inp.write("ATOMIC_POSITIONS ")
      if self.content['system']['ibrav'] == 0:
        inp.write("angstrom\n")
        R = molecule.R
      else:
        inp.write("\n")
        R = molecule.R_scale
      for a in range(len(type_list)):
        inp.write(' %-3s' % type_list[a])
        for i in range(3):
          inp.write(' % 12.8f' % R[a, i])
        inp.write("\n")
      inp.write("\n")

      if 'kmesh' in setting and setting['kmesh']:
        inp.write("K_POINTS automatic\n")
        for k in setting['kmesh']:
          inp.write(" %d" % k)
        for s in range(3):
          inp.write(" 0")
        inp.write('\n\n')

      if 'save_restart' in setting and setting['save_restart']:
        inp.write("ALCHEMY reference\n\n")

      if 'restart' in setting and setting['restart']:
        if 'scf_step' in setting and setting['scf_step'] == 1:
          inp.write("ALCHEMY prediction\n\n")

      if self.content['system']['ibrav'] == 0:
        inp.write("CELL_PARAMETERS angstrom\n")
        if 'lattice' not in self.setting:
          self.celldm2lattice()
        lattice_vec = self.setting['lattice']
        for vec in lattice_vec:
          for component in vec:
            inp.write(' %9.6f' % component)
          inp.write('\n')
  
      for pp in pp_files:
        pp_file = os.path.join(qtk.setting.espresso_pp, pp)
        if pp not in inp.dependent_files:
          inp.dependent_files.append(pp_file)

      if 'no_cleanup' in setting and setting['no_cleanup']:
        inp.close(no_cleanup=True)
      else:
        inp.close()

      return inp

    if 'ref_dir' in self.setting:
      inp = refCopy(name, **self.setting)
    else:
      setting = copy.deepcopy(self.setting)
      inp = writeInp(name, **setting)

    return inp

def PPString(inp, mol, i, n, outFile):
  """
  append PP file names to inp.pp_files
  """
  alchemy = re.compile('^\w*2\w*_\d\d\d$')
  ppstr = re.sub('\*', '', mol.string[i])
  if ppstr:
    PPStr = ppstr
    pp_root, pp_ext = os.path.split(ppstr)
  else:
    if inp.setting['pp_type'] == 'geodecker':
      element = mol.type_list[i].title()
      if 'd_shell' in inp.setting:
        if type(inp.setting['d_shell']) is not list:
          inp.setting['d_shell'] = [inp.setting['d_shell']]
      if qtk.n2ve(mol.type_list[i].title()) > 10:
        shell = '-d'
      elif 'd_shell' in inp.setting \
      and element in inp.setting['d_shell']:
        shell = '-d'
      else:
        shell = ''
      pp_xc_dict = {
        'lda': 'pz',
        'pbe0': 'pbe',
        'b3lyp': 'blyp',
      }
      pp_xc = inp.setting['pp_theory'].lower()
      if pp_xc in pp_xc_dict:
        pp_xc = pp_xc_dict[pp_xc]
      PPStr = ''.join([c for c in mol.type_list[i] if not c.isdigit()])\
              + '.' + pp_xc + shell + '-hgh.UPF'
    elif inp.setting['pp_type'] == 'cpmd':
      PPStr = PPName(inp, mol, i, n)
  xc = inp.setting['pp_theory'].lower()
  if not mol.string[i]:
    if inp.setting['pp_type'] == 'geodecker':
      PPCheck(pp_xc, mol.type_list[i].title(), PPStr)
    elif inp.setting['pp_type'] == 'cpmd':
      saved_pp = PPCheck_cpmd(pp_xc, mol.type_list[i].title(), PPStr)
      new_pp1 = saved_pp + '.UPF'
      conv_pp = sp.Popen("%s %s" % \
        (qtk.setting.espresso_cpmd2upf_exe, saved_pp),
        shell=True)
      conv_pp.wait()
      new_pp1_file = os.path.split(new_pp1)[1]
      new_pp1_trg = os.path.join(qtk.setting.espresso_pp, new_pp1_file)
      if not os.path.exists(new_pp1_trg):
        shutil.copy(new_pp1, qtk.setting.espresso_pp)
      PPStr = PPStr + '.UPF'

  elif alchemy.match(mol.string[i]):
    cpmd_pp = alchemyPP(xc, PPStr)
    new_pp1 = cpmd_pp + '.UPF'
    if not os.path.exists(new_pp1):
      qtk.report('espresso', "rewrite Goedecker's PP to UPF")
      conv_pp = sp.Popen("%s %s" % \
        (qtk.setting.espresso_cpmd2upf_exe, cpmd_pp),
        shell=True)
      conv_pp.wait()
      if conv_pp.returncode != 0:
        # dirty fix for espresso alchemy conversion routine
        qtk.warning('conversion failed..., trying path end points')
        root, _ = os.path.splitext(PPStr)
        element_str = re.sub('_.*', '', root)
        element1 = re.sub('2.*', '', element_str)
        element2 = re.sub('.*2', '', element_str)
        fraction = float(re.sub('.*_', '', root))/100
        if fraction == 0.0:
          strpp = element1 + "_q" + str(qtk.n2ve(element1)) +\
                  "_" + xc + '.psp'
        elif fraction == 1.0:
          strpp = element2 + "_q" + str(qtk.n2ve(element2)) +\
                  "_" + xc + '.psp'
        else:
          qtk.exit("PP conversion failed for intermediate lambda")
        strpp = os.path.join(qtk.setting.cpmd_pp, strpp)
        conv_pp = sp.Popen("%s %s" % \
          (qtk.setting.espresso_cpmd2upf_exe, strpp),
          shell=True)
        conv_pp.wait()
        os.rename(strpp + '.UPF', new_pp1)
    new_pp1_file = os.path.split(new_pp1)[1]
    new_pp1_trg = os.path.join(qtk.setting.espresso_pp, new_pp1_file)
    if not os.path.exists(new_pp1_trg):
      shutil.copy(new_pp1, qtk.setting.espresso_pp)
    PPStr = PPStr + '.UPF'

  return PPStr

def PPCheck(xc, element, pp_file_str, **kwargs):
  ne = qtk.n2ve(element)
  saved_pp_path = os.path.join(qtk.setting.espresso_pp, pp_file_str)
  if not os.path.exists(saved_pp_path) and qtk.setting.download_pp:
    url = os.path.join(qtk.setting.espresso_pp_url, pp_file_str)
    try:
      pp_content = urllib2.urlopen(url).read()
      qtk.report('', 'pp file %s not found in %s. ' \
                 % (pp_file_str, qtk.setting.espresso_pp) + \
                 'but found in espresso page, download now...')
      new_pp_file = open(saved_pp_path, 'w')
      new_pp_file.write(pp_content)
      new_pp_file.close()
    except:
      qtk.warning('something wrong with pseudopotential')

class out(PlanewaveOutput):
  def __init__(self, qmout, **kwargs):
    PlanewaveOutput.__init__(self, qmout, **kwargs)
    out_file = open(qmout)
    data = out_file.readlines()
    out_file.close()
    Et_pattern = re.compile("^.*total energy *=.*$")
    Et_list = filter(Et_pattern.match, data)
    if len(Et_list) > 0:
      Et_str = filter(Et_pattern.match, data)[-1]
      Et = float(Et_str.split()[-2])
      self.Et, self.unit = qtk.convE(Et, 'Ry-Eh')
      out_folder = os.path.split(os.path.abspath(qmout))[0]
      save = glob.glob(os.path.join(out_folder, '*.save'))
      # extract band structure from xml files
      if save:
        save = save[0]
        try:
          data_xml = os.path.join(save, 'data-file.xml')
          xml_file = open(data_xml)
          tree = ET.parse(xml_file)
          xml_file.close()
          self.xml = tree.getroot()
          kpoints = []
          band = []
          # access data for each kpoint
          for k in self.xml[-2]:
            k_str = k[0].text
            coord = [float(c) for c in k_str.split()]
            weight = float(k[1].text.split()[0])
            coord.append(weight)
            kpoints.append(coord)
            ev_file = os.path.join(save, k[2].attrib['iotk_link'])
            k_xml_file = open(ev_file)
            k_xml = ET.parse(k_xml_file)
            k_xml_file.close()
            ev_k = k_xml.getroot()
            ev_str = ev_k[2].text.split()
            ev = [qtk.convE(float(entry), 'Eh-eV')[0]\
                  for entry in ev_str]
            band.append(ev)
            occ_str = ev_k[3].text.split()
            occ = [float(entry) for entry in occ_str]
          self.kpoints = np.array(kpoints)
          self.mo_eigenvalues = copy.deepcopy(band[0])
          self.band = np.array(band)
          self.occupation = occ
          diff = np.diff(occ)
          pos = diff[np.where(abs(diff) > 0.5)]
          mask = np.in1d(diff, pos)
          ind = np.array(range(len(diff)))
          if len(ind[mask]) > 0:
            N_state = ind[mask][0]
            vb = max(self.band[:, N_state])
            cb = min(self.band[:, N_state + 1])
            vb_pos = np.argmax(self.band[:, N_state])
            cb_pos = np.argmin(self.band[:, N_state + 1])
            self.Eg = cb - vb
            if vb_pos == cb_pos:
              self.Eg_direct = True
            else:
              self.Eg_direct = False
            
        except IOError:
          qtk.warning('xml file of job %s not found' % qmout)
    else:
      qtk.warning('job %s not finished' % qmout)
