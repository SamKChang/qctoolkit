import qctoolkit, re, sys
import setting
import numpy as np
from qctoolkit import utilities as ut

class inp(setting.QMSetting):
  def __init__(self, structure_inp, info):
    self.setting = setting.QMSetting()
    self.set_center = False
    self.set_celldm = False
    self.set_margin = False
    self.set_mode = False
    self.set_step = False
    self.debug = False
    self.restart = False

    self.atom_list = {}

    if structure_inp.endswith("xyz"):
      self.name = re.sub(".xyz", "", structure_inp)
      self.structure = qctoolkit.geometry.Molecule()
      self.structure.read_xyz(structure_inp)
    else:
      print "unknown structure format"
    self.info = info
   
  def PPStringDefault(self, atom_type):
    return atom_type.title() + "_q" + str(ut.n2ve(atom_type))\
           + "_" + self.setting.theory.lower() + ".psp"

  # CPMD input format
  def write(self, name):
    inp= sys.stdout if re.match("stdout", name) else open(name,"w")

    isolated = False
    if re.match("ISOLATED", self.setting.symmetry.upper()):
      isolated = True
  
    angstrom = False
    if re.match("ANGSTROM", self.setting.unit.upper()):
      angstrom = True
  
    print >>inp, "&INFO"
    print >>inp, " qmInfo:%s" % self.info
    print >>inp, "&END"
    print >>inp, ""
    print >>inp, "&CPMD"

    if self.debug:
      print >>inp, " BENCHMARK"
      print >>inp, "  1 0 0 0 0 0 0 0 0 0 "

    if not self.set_mode:
      print >>inp, " OPTIMIZE WAVEFUNCTION"
    elif re.match("GEOPT", self.setting.mode.upper()):
      print >>inp, " OPTIMIZE GEOMETRY"
    elif re.match("KSEG", self.setting.mode.upper()):
      print >>inp, " KOHN-SHEM ENERGIES"
      print >>inp, "  " + str(self.setting.ks_states)
      self.restart = True
    else:
      sys.exit("ERROR from io_format/cpmd.py->inp.write: " +\
               "mode '" + self.setting.mode +\
               "' is not implemented. Supported modes:\n" +\
               " single_point\n" +\
               " geopt\n" +\
               " kseg"
              )

    if self.restart:
      print >>inp, " RESTART WAVEFUNCTION"

    if not self.set_center and self.set_margin:
      print >>inp, " CENTER MOLECULE OFF"
      if not self.set_celldm:
        Rt = np.transpose(self.structure.R)
        for i in xrange(0, 3):
          self.setting.celldm[i] = (max(Rt[i]) - min(Rt[i]))\
                                 + 2 * self.setting.margin
        new_center=[min(Rt[i])-self.setting.margin \
          for i in (0,1,2)]
        self.structure.center(new_center)
      else:
        sys.exit("ERROR from io_format/cpmd.py->inp.write: " +\
                 "celldm and margin " + \
                 "can NOT be set simultaneously.")
  
    elif self.set_center and not self.set_margin:
      self.structure.center(self.setting.center)
      print >>inp, " CENTER MOLECULE OFF"

    elif self.set_center and self.set_margin:
      sys.exit("ERROR from io_format/cpmd.py->inp.write: " +\
               "center and margin " + \
               "can NOT be set simultaneously.")

    if self.set_step:  
      print >>inp, " MAXITER"
      print >>inp, "  " + str(self.setting.maxstep)

    if self.setting.save_density:
      print >>inp, " RHOOUT"

    if self.structure.multiplicity > 1:
      print >>inp, " LOCAL SPIN DENSITY"
    print >>inp, " MIRROR"
    print >>inp, "&END"
    print >>inp, ""

    print >>inp, "&DFT"
    print >>inp, " FUNCTIONAL " + self.setting.theory.upper()
    print >>inp, "&END"
        
    print >>inp, ""
    print >>inp, "&SYSTEM"
    print >>inp, " SYMMETRY"
    print >>inp, "  " + self.setting.symmetry.upper()
    if isolated:
      print >>inp, " POISSON SOLVER TUCKERMAN"
    if angstrom:
      print >>inp, " ANGSTROM"
    print >>inp, " CELL ABSOLUTE"
    print >>inp, "  " + " ".join(map(str, self.setting.celldm))
    print >>inp, " CUTOFF"
    print >>inp, "  " + str(self.setting.cutoff)
    if self.setting.mesh:
      mesh = np.array(self.setting.celldm[0:3])\
           * self.setting.mesh
      mesh.tolist()
      print >>inp, " MESH"
      print >>inp, "  " + " ".join(map(str,mesh))
    if not isolated:
      print >>inp, " KPOINTS MONKHORST-PACK"
      print >>inp, "  " + " ".join(map(str, self.setting.kmesh))
    if self.structure.charge:
      print >>inp, " CHARGE"
      print >>inp, "  " + str(self.structure.charge)
    if self.structure.multiplicity > 1:
      print >>inp, " MULTIPLICITY"
      print >>inp, "  " + str(self.structure.multiplicity)
    print >>inp, "&END"
    print >>inp, ""
    print >>inp, "&ATOMS"
  
    # loop through all atom types
    self.structure.sort()
    PP = self.PPStringDefault
    type_index = self.structure.index
    type_list = self.structure.type_list
    atom_list = self.atom_list
    for atom_type in xrange(0,len(type_index)-1):
      type_n = type_index[atom_type+1] - type_index[atom_type]
      if atom_list.has_key(type_list[type_index[atom_type]]):
        key = type_list[type_index[atom_type]]
        print >>inp, "*" + atom_list[key]
      else:
        print >>inp, "*" + \
                PP(self.structure.type_list[type_index[atom_type]])
      print >>inp, " LMAX=F\n" + "  " + str(type_n)
      for I in\
        xrange(type_index[atom_type],type_index[atom_type+1]):
        print >>inp, "  " + \
        " ".join(" % 8.4f" % x for x in self.structure.R[I][:])
      print >>inp 
    print >>inp, "&END"

    if not re.match("stdout", name):
      inp.close()

class out(object):
  def __init__(self, qmout):
    self.info = ''
    self.getEt(qmout)
    #self.getSteps(qmout)

  def getEt(self, name):
    out = sys.stdout if re.match('stdout',name)\
          else open(name, 'r')

    done = False
    finished = False
  
    scf_p = re.compile('^ *[0-9]*  [0-9]\.[0-9]{3}E-[0-9]{2}   .*')
    Et_cpmd = re.compile('.*TOTAL ENERGY = *([-0-9\.]*)')
    done_cpmd = re.compile(' \* *FINAL RESULTS *\*')
    qmInfo = re.compile('.*qmInfo.*')
    info_head = re.compile('.*qmInfo:')
    info_tail = re.compile(' *\*\*.*')

    while True:
      line = out.readline()
      if not line: break

      if (re.match(scf_p, line)):
        data = [float(x) for x in line.split()]
        self.SCFStep = int(data[0])
      elif re.match(qmInfo, line):
        tmp1 = re.sub(info_head, '', line)
        tmp2 = re.sub(info_tail, '', tmp1)
        self.info = tmp2
      elif (re.match(done_cpmd,line)):
        done = True,
      elif (re.match(Et_cpmd, line)) and done:
        self.Et = float(Et_cpmd.match(line).group(1))
        finished = True
        
    if not finished:
      self.Ehartree = 0
      self.status = 'not finished'
    out.close()


  def getSteps(self, name):
    # CPMD format
    scf_p = \
      re.compile('^ *[0-9]*  [0-9]\.[0-9]{3}E-[0-9]{2}   .*')
    out = open(name, 'r')
    while True:
      line = out.readline()
      if not line: break
      if (re.match(scf_p, line)):
        data = [float(x) for x in line.split()]
        self.SCFStep = int(data[0])
    #return SCFStep
