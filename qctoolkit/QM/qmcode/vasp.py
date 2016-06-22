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
    PlanewaveInput.__init__(self, molecule, **kwargs)
    self.setting['pp_type'] = 'vasp_default'
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
      self.pp_path = qtk.setting.vasp_pp
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
    if self.setting['theory'] == 'pbe0':
      incar.write("LHFCALC = .TRUE.\n")
      incar.write("GGA = PE\n")
    elif self.setting['theory'] == 'hse06':
      incar.write("LHFCALC = .TRUE.\n")
      incar.write("HFSCREEN = 0.2 \n")
      incar.write("GGA = PE\n")

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
      xml_file = open(qmoutXML)
      tree = ET.parse(qmoutXML)
      xml_file.close()
      self.xml = tree.getroot()
      self.Et, self.unit = qtk.convE(float(self.xml[-2][-5][1].text),
                                     'eV-Eh', '-')
      self.info = self.xml[0][1].text
      self.SCFStep = len(self.xml[-2])-9
