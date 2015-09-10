import qctoolkit as qtk
import qctoolkit.io_format.setting_pw as pw
import qctoolkit.io_format.qminp as qin
import os, sys, copy

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

    def catPOTCAR(PPPath):
      print PPPath

    def n_list(atom_number):
      print atom_number

    def R_list(coord):
      print " ".join(" % 8.4f" % x for x in coord)


    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # write to INCAR and generate POTCAR
    qtk.report("vasp.inp", "writing", path+"INCAR")
    new_structure.sort()
    type_index = new_structure.index
    type_list = new_structure.type_list
    atom_list = self.atom_list
    Z = new_structure.Z
    for atom_type in xrange(0,len(type_index)-1):
      type_n = type_index[atom_type+1] - type_index[atom_type]
      if atom_list.has_key(str(Z[type_index[atom_type]])):
        key = str(Z[type_index[atom_type]])
        AtomPP = os.path.join(self.PP, atom_list[key], 'POTCAR')
        catPOTCAR(AtomPP)
        del atom_list[str(Z[type_index[atom_type]])]
      else:
        type_name = qtk.Z2n(Z[type_index[atom_type]])
        AtomPP = os.path.join(self.PP, type_name, 'POTCAR')
        catPOTCAR(AtomPP)
      n_list(type_n)
      for I in\
        xrange(type_index[atom_type],type_index[atom_type+1]):
        R_list(new_structure.R[I][:])
      print >>incar





    # !!!!!!!!!!!!!!!!
    # write to KPOINTS
    qtk.report("vasp.inp", "writin", path+"KPOINTS")
    print >> kpoints, "kp yo"

    # !!!!!!!!!!!!!!!
    # write to POSCAR
    qtk.report("vasp.inp", "writing", path+"POSCAR")
    print >> poscar, "poscar yo"

    if name:
      os.chdir(cwd)
      incar.close()
      kpoints.close()
      poscar.close()

class out(object):
  def __init__(self, qmout):
    pass
