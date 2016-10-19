import qctoolkit as qtk
from qctoolkit.QM.planewave_io import PlanewaveInput
from qctoolkit.QM.planewave_io import PlanewaveOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import xml.etree.ElementTree as ET
import universal as univ

class inp(PlanewaveInput):
  def __init__(self, molecule, **kwargs):
    if 'cutoff' not in kwargs:
      cutoff_default = True
    else:
      cutoff_default = False
    PlanewaveInput.__init__(self, molecule, **kwargs)
    if cutoff_default:
      del self.setting['cutoff']
    self.setting['pp_type'] = 'PAW'
    self.setting.update(kwargs)
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
    self.setting['root_dir'] = name
    self.setting['no_molecule'] = False
    if name:
      self.setting['output'] = True
      name = os.path.splitext(name)[0]
    else: 
      self.setting['output'] = False
    incar, molecule = \
      super(PlanewaveInput, self).write('INCAR', **self.setting)
    self.setting['no_molecule'] = True
    kpoints = \
      super(PlanewaveInput, self).write('KPOINTS', **self.setting)
    poscar = \
      super(PlanewaveInput, self).write('POSCAR', **self.setting)
    potcar = \
      super(PlanewaveInput, self).write('POTCAR', **self.setting)

    # !!!!!!!!!!! TODO !!!!!!!!!!!
    # Center molecule
    # charge multiplicity
    # optimizer
    # 
    # write CPMD to modulize structure manipulation?

    PPPath = []
    n_list = []
    R_list = []

    def catPOTCAR(path):
      if os.path.exists(path):
        PPPath.append(path)
      else:
        qtk.exit("PP file: " + path + " not found")

    def getNlist(atom_number):
      n_list.append(atom_number)

    def getRlist(coord):
      R_list.append(coord)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # construct to POSCAR and POTCAR data
    molecule.sort()
    type_index = molecule.index
    type_list = molecule.type_list
    Z = molecule.Z

    self.pp_path = None
    if 'pp_path' not in self.setting:
      if 'pp_theory' not in self.setting:
        theory_dict = {
          'pbe': 'pbe',
          'pbe0': 'pbe',
          'hse06': 'pbe',
          'hse03': 'pbe',
          'lda': 'lda',
        }
        if self.setting['theory'] not in theory_dict:
          qtk.warning("%s is not supported, change theory to LDA" \
                      % (self.setting['theory']))
          self.setting['theory'] = 'lda'
        theory = theory_dict[self.setting['theory']]
        if theory.lower() not in ['pbe', 'lda']:
          qtk.warning('xc: %s is not supported, using LDA PP' % \
                      theory.upper())
          theory = 'LDA'
        self.pp_path = qtk.setting.vasp_pp + '_%s_%s' % \
                       (theory.upper(), 
                        self.setting['pp_type'].upper())
      else:
        self.pp_path = qtk.setting.vasp_pp + '_%s_%s' % \
                       (self.setting['pp_theory'].upper(), 
                        self.setting['pp_type'].upper())
    else:
      self.pp_path = self.setting['pp_path']

    for atom_type in xrange(0,len(type_index)-1):
      type_n = type_index[atom_type+1] - type_index[atom_type]
      # check special PP folder
        # not yet implemented
      # default PP path
      type_name = type_list[type_index[atom_type]]
      AtomPP = os.path.join(self.pp_path, type_name, 'POTCAR')
      catPOTCAR(AtomPP)
      getNlist(type_n)
      for I in\
        xrange(type_index[atom_type],type_index[atom_type+1]):
        getRlist(molecule.R[I][:])

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # write to INCAR and generate POTCAR
    qtk.report("vasp.inp", "writing", "INCAR")
    incar.write("SYSTEM = %s\n" % self.setting['info'])
    incar.write("ISMEAR = 0\n")
    incar.write("IBRION = 2\n")
    if 'restart' in self.setting and self.setting['restart']:
      incar.write("ISTART = 1\n")
    if 'cutoff' in self.setting:
      cutoff = self.setting['cutoff']
      incar.write("ENCUT = %.2f" % (cutoff * 13.605698066))
      incar.write(" # in eV, that is %.1f Ry\n" % cutoff)
    if 'scf_step' in self.setting:
      incar.write('NELM = %d\n' % self.setting['scf_step'])
    if 'vdw' in self.setting:
      vdw = self.setting['vdw'].lower()
      if vdw != 'none':
        if vdw=='d2':
          incar.write("IVDW = 10\n")
        elif vdw=='d3':
          incar.write("IVDW = 11\n")
        elif vdw=='d3-bj':
          incar.write("IVDW = 12\n")
        elif vdw=='mbd':
          incar.write("IVDW = 202\n")
        elif vdw=='mbd_iter':
          incar.write("IVDW = 212\n")
        else:
          qtk.exit("VDW '%s' is not supported for VASP" % vdw)
    if 'ks_states' in self.setting:
      vs = int(round(self.molecule.getValenceElectrons() / 2.0))
      nbnd = vs + self.setting['ks_states']
      incar.write("NBANDS = %d\n" % nbnd)
    if 'full_kmesh' in self.setting and self.setting['full_kmesh']:
      incar.write("ISYM = -1\n")
    if self.setting['theory'] == 'pbe0':
      incar.write("LHFCALC = .TRUE.\n")
      incar.write("GGA = PE\n")
    elif self.setting['theory'] == 'hse06':
      incar.write("GGA = PE\n")
      incar.write("\n##HSE setting\n")
      incar.write("LHFCALC = .TRUE.\n")
      incar.write("HFSCREEN = 0.2 \n")
      incar.write("ALGO = D\n")

    if molecule.charge != 0:
      nve = molecule.getValenceElectrons()
      incar.write("NELECT = %d\n" % (nve))
    if 'save_density'not in self.setting\
    or not self.setting['save_density']:
      incar.write("LCHARG = .FALSE.\n")
    if 'scalapack' not in self.setting:
      incar.write("LSCALAPACK = .FALSE.\n")
    elif not self.setting['scalapack']:
      incar.write("LSCALAPACK = .FALSE.\n")
    incar.close()

    # !!!!!!!!!!!!!!!!
    # write to KPOINTS
    qtk.report("vasp.inp", "writing", "KPOINTS")
    if 'kmesh' not in self.setting:
      kpoints.write("Gamma-point only\n")
      kpoints.write(" 1       ! one k-point\n")
      kpoints.write("rec      ! in units of reciprocal vector\n")
      kpoints.write(" 0 0 0 1 ! coordinates and weight\n")
    else:
      k1, k2, k3 = self.setting['kmesh']
      kpoints.write("Automatic mesh\n")
      kpoints.write(" 0       ! number of k-points = 0")
      kpoints.write(" ->automatic generation scheme\n")
      kpoints.write("Gamma    ! generate a Gamma centered grid\n")
      kpoints.write(" %d %d %d   ! number of k grids\n" % (k1, k2, k3))
    kpoints.close(no_cleanup=True)

    # !!!!!!!!!!!!!!!
    # write to POSCAR
    qtk.report("vasp.inp", "writing", "POSCAR")
    poscar.write(self.setting['info'] + '\n')
    poscar.write("1.0\n")
    self.celldm2lattice()
    for i in range(3):
      for j in range(3):
        poscar.write(" %7.4f" % self.setting['lattice'][i,j])
      poscar.write(" ! lattic vector a(%d)\n" %i)
    for n in n_list:
      poscar.write(str(n) + ' ')
    poscar.write("! number of atoms in order of POTCAR\n")
    poscar.write("cart ! cartesian coordinates\n")
    for R in R_list:
      for X in R:
        poscar.write(" %7.4f" % X)
      poscar.write("\n")
    poscar.close(no_cleanup=True)

    # !!!!!!!!!!!!!!!
    # write to POTCAR
    qtk.report("vasp.inp", "writing", "POTCAR")
    for PP_file in PPPath:
      qtk.report("vasp.inp.POTCAR", PP_file)
      if name:
        with open(PP_file) as PP:
          for line in PP:
            potcar.write(str(line))
      else:
        potcar.write("cat %s\n" % PP_file)
    potcar.close(no_cleanup=True)

    return incar

class out(PlanewaveOutput):
  """
  directly parse vasp xml output, 'vasprun.xml'
  converged energy, system info, and scf steps are extracted
  """
  def __init__(self, qmoutXML, **kwargs):
    PlanewaveOutput.__init__(self, qmoutXML, **kwargs)
    if qmoutXML:
      #xml_file = open(qmoutXML)
      tree = ET.parse(qmoutXML)
      #xml_file.close()
      self.xml = tree.getroot()
      self.Et, self.unit = qtk.convE(float(self.xml[-2][-5][1].text),
                                     'eV-Eh', '-')
      self.info = self.xml[0][1].text
      self.SCFStep = len(self.xml[-2])-9

      try:
        # kpoints and band structure data
        kpoints = []
        band = []
        occ = []
  
        for i in range(len(self.xml[2][1])):
          k_str = self.xml[2][1][i].text
          weight = float(self.xml[2][2][i].text)
          band_k = []
          occ_k = []
          bk_xml = self.xml[-2][-3][0][5][0][i]
          for b_xml in bk_xml:
            b, o = [float(c) for c in b_xml.text.split()]
            band_k.append(b)
            occ_k.append(o)
          coord = [float(c) for c in k_str.split()]
          coord.append(weight)
          kpoints.append(coord)
          band.append(band_k)
          occ.append(occ_k)
        self.mo_eigenvalues = copy.deepcopy(band[0])
        self.kpoints = np.array(kpoints)
        self.band = np.array(band)
        self.occupation = occ[0]
  
        diff = np.diff(occ[0])
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
  
        # DOS data
        dos = []
        self.E_fermi = float(self.xml[-2][-1][0].text) 
        print self.E_fermi
        for dos_xml in self.xml[-2][-1][1][0][5][0]:
          dos_lst = filter(None, dos_xml.text.split(' '))
          dos_E = [float(d) for d in dos_lst]
          dos.append(dos_E)
        self.dos = np.array(dos)

      except:
        qtk.warning("error when accessing kpoint data for %s"\
                    % self.name)
