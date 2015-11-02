import qctoolkit as qtk
from qctoolkit.QM.planewave_io import PlanewaveInput
from qctoolkit.QM.planewave_io import PlanewaveOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import xml.etree.ElementTree as ET

class inp(PlanewaveInput):
  def __init__(self, molecule, **kwargs):
    PlanewaveInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)

  def run(self, name=None, **kwargs):
    if not name:
      name = self.molecule.name
    cwd = os.getcwd()
    self.write(name)
    os.chdir(name)
    try:
      out = qmjob.QMRun(name, 'vasp', **kwargs)
    except:
      qtk.warning("qmjob finished unexpectedly for '" + \
                  name + "'")
      out = None
    finally:
      os.chdir(cwd)
    return out
    
  def write(self, name=None):
    molecule = copy.deepcopy(self.molecule)
    self.cm_check(molecule)
    if name:
      cwd = os.getcwd()
      if os.path.exists(name) and not qtk.setting.no_warning:
        qtk.prompt(name + ' exists, overwrite?')
        try:
          shutil.rmtree(name)
        except:
          qtk.exit("can not remove folder: " + name)
      os.mkdir(name)
      os.chdir(name)

    incar   = open('INCAR',  'w') if name else sys.stdout
    kpoints = open('KPOINTS','w') if name else sys.stdout
    poscar  = open('POSCAR', 'w') if name else sys.stdout
    potcar  = open('POTCAR', 'w') if name else sys.stdout

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
        qtk.exit("PP file: " + path + "not found")

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

    if not qtk.setting.vasp_pp:
      qtk.exit("vasp pp_path not set")
    else:
      self.pp_path = qtk.setting.vasp_pp

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
    print >> incar, "SYSTEM = %s" % self.setting['info']
    print >> incar, "ISMEAR = 0"
    print >> incar, "IBRION = 2"
    if 'scf_step' in self.setting:
      incar.write('NELM = %d\n' % self.setting['scf_step'])
    if 'vdw' in self.setting:
      vdw = self.setting['vdw'].lower()
      if vdw != 'none':
        if vdw=='d2':
          print >> incar, "IVDW = 10"
        elif vdw=='d3':
          print >> incar, "IVDW = 11"
        elif vdw=='d3-bj':
          print >> incar, "IVDW = 12"
        elif vdw=='mbd':
          print >> incar, "IVDW = 202"
        elif vdw=='mbd_iter':
          print >> incar, "IVDW = 212"
        else:
          qtk.exit("VDW '%s' is not supported for VASP" % vdw)
    if molecule.charge != 0:
      nve = molecule.getValenceElectrons()
      print >> incar, "NELECT = %d" % (nve)
    if 'save_density'not in self.setting\
    or self.setting['save_density']:
      print >> incar, "LCHARG = .FALSE."
    if 'scalapack' not in self.setting:
      print >> incar, "LSCALAPACK = .FALSE."
    elif not self.setting['scalapack']:
      print >> incar, "LSCALAPACK = .FALSE."

    # !!!!!!!!!!!!!!!!
    # write to KPOINTS
    qtk.report("vasp.inp", "writin", "KPOINTS")
    if 'kmesh' not in self.setting:
      print >> kpoints, "Gamma-point only"
      print >> kpoints, " 1       ! one k-point"
      print >> kpoints, "rec      ! in units of reciprocal vector"
      print >> kpoints, " 0 0 0 1 ! coordinates and weight"

    # !!!!!!!!!!!!!!!
    # write to POSCAR
    qtk.report("vasp.inp", "writing", "POSCAR")
    print >> poscar, self.setting['info']
    print >> poscar, "1.0"
    self.celldm2lattice()
    for i in range(3):
      for j in range(3):
        print >> poscar, " %7.4f" % self.setting['lattice'][i,j],
      print >> poscar, "! lattic vector a(%d)" %i
    for n in n_list:
      print >> poscar, n,
    print >> poscar, "! number of atoms in order of POTCAR"
    print >> poscar, "cart ! cartesian coordinates"
    for R in R_list:
      for X in R:
        print >> poscar, " %7.4f" % X,
      print >> poscar

    # !!!!!!!!!!!!!!!
    # write to POTCAR
    qtk.report("vasp.inp", "writing", "POTCAR")
    for PP_file in PPPath:
      qtk.report("vasp.inp.POTCAR: ", PP_file)
      if name:
        with open(PP_file) as PP:
          for line in PP:
            print >> potcar, line,
      else:
        print >> potcar, "cat %s" % PP_file

    if name:
      incar.close()
      kpoints.close()
      poscar.close()
      potcar.close()

    os.chdir(cwd)

class out(PlanewaveOutput):
  """
  directly parse vasp xml output, 'vasprun.xml'
  converged energy, system info, and scf steps are extracted
  """
  def __init__(self, qmoutXML, **kwargs):
    PlanewaveOutput.__init__(self, qmoutXML, **kwargs)
    tree = ET.parse(qmoutXML)
    self.xml = tree.getroot()
    self.Et = qtk.convE(float(self.xml[-2][-5][1].text), 'eV-Eh')
    self.info = self.xml[0][1].text
    self.SCFStep = len(self.xml[-2])-9
