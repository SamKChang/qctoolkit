import qctoolkit as qtk
import qctoolkit.io_format.setting_pw as pw
import qctoolkit.io_format.qminp as qin
import os, sys, copy, shutil

def qmDir():
  pass

class inp(qin.PwInp):
  def __init__(self, structure_inp, info, **kwargs):
    qin.PwInp.__init__(self, structure_inp, info, **kwargs)
    if 'PP' in kwargs:
      self.PP = kwargs['PP']
    else:
      self.PP = qtk.PP

  def write(self, *args, **kwargs):
    new_structure = copy.deepcopy(self.structure)
    if len(args) == 1:
      name = args[0]
      if os.path.exists(name):
        if 'new_run' in kwargs and kwargs['new_run']:
          shutil.rmtree(name)
          cwd = os.getcwd()
          os.mkdir(name)
          os.chdir(name)
          path = name + '/'
        else:
          qtk.warning("inp.write: path " + name + " exist, "+\
                      "nothing to be done")
      else:
        cwd = os.getcwd()
        os.mkdir(name)
        os.chdir(name)
        path = name + '/'
    else: 
      name = ''
      path = ''
    incar   = open('INCAR', 'w')   if name else sys.stdout
    kpoints = open('KPOINTS', 'w') if name else sys.stdout
    poscar  = open('POSCAR', 'w')  if name else sys.stdout
    potcar  = open('POTCAR', 'w')  if name else sys.stdout


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
#      print " ".join(" % 8.4f" % x for x in coord)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # construct to POSCAR and POTCAR data
    new_structure.sort()
    type_index = new_structure.index
    type_list = new_structure.type_list
    atom_list = self.atom_list
    Z = new_structure.Z
    for atom_type in xrange(0,len(type_index)-1):
      type_n = type_index[atom_type+1] - type_index[atom_type]
      # check special PP folder
      if atom_list.has_key(str(Z[type_index[atom_type]])):
        key = str(Z[type_index[atom_type]])
        AtomPP = os.path.join(self.PP, atom_list[key], 'POTCAR')
        catPOTCAR(AtomPP)
        del atom_list[str(Z[type_index[atom_type]])]
      # default PP path
      else:
        type_name = qtk.Z2n(Z[type_index[atom_type]])
        AtomPP = os.path.join(self.PP, type_name, 'POTCAR')
        catPOTCAR(AtomPP)
      getNlist(type_n)
      for I in\
        xrange(type_index[atom_type],type_index[atom_type+1]):
        getRlist(new_structure.R[I][:])

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # write to INCAR and generate POTCAR
    qtk.report("vasp.inp", "writing", path+"INCAR")
    print >> incar, "SYSTEM = %s" % self.info
    print >> incar, "ISMEAR = 0"
    print >> incar, "IBRION = 2"

    # !!!!!!!!!!!!!!!!
    # write to KPOINTS
    qtk.report("vasp.inp", "writin", path+"KPOINTS")
    if not self.setting.kpoints:
      print >> kpoints, "Gamma-point only"
      print >> kpoints, " 1       ! one k-point"
      print >> kpoints, "rec      ! in units of reciprocal vector"
      print >> kpoints, " 0 0 0 1 ! coordinates and weight"

    # !!!!!!!!!!!!!!!
    # write to POSCAR
    qtk.report("vasp.inp", "writing", path+"POSCAR")
    print >> poscar, self.info
    print >> poscar, "1.0"
    for i in range(3):
      for j in range(3):
        print >> poscar, " %7.4f" % self.setting.lattice[i,j],
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
    qtk.report("vasp.inp", "writing", path+"POTCAR")
    for PP_file in PPPath:
      qtk.report("vasp.inp.POTCAR: ", PP_file)
      if name:
        with open(PP_file) as PP:
          for line in PP:
            print >> potcar, line,
      else:
        print >> potcar, "cat %s" % PP_file
  

    if name:
      os.chdir(cwd)
      incar.close()
      kpoints.close()
      poscar.close()
      potcar.close()

class out(object):
  def __init__(self, qmout):
    pass
