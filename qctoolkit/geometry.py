import numpy as np
import utilities as ut
import re
import sys
import copy

class Molecule(object):
  def __init__(self):
    # number of atoms
    self.N = 0
    # atom coordinates
    self.R = np.atleast_2d(np.array([]))
    # atom symbols
    self.type_list = 'None'
    # nuclear charges
    self.Z = 0
    # moelcule charge
    self.charge = 0
    self.multiplicity = 1
    # index of different atoms
    self.index = 0
    self.bonds = {}
    self.bond_types = {}

#  def have_bond(self, type_a, type_b):
#    obmol = ut.qt2ob(self)

  def __add__(self, other):
    out = Molecule()
    out.N = self.N + other.N
    out.R = np.vstack([self.R, other.R])
    out.Z = np.hstack([self.Z, other.Z])
    return out

  def find_bonds(self):
    itr = 0
    for i in xrange(self.N):
      for j in xrange(i+1, self.N):
        d_ij = np.linalg.norm(self.R[i,:] - self.R[j,:])
        if d_ij < 1.75:
          if self.Z[i] < self.Z[j]:
            atom_begin = self.Z[i]
            atom_end = self.Z[j]
            index_begin = i
            index_end = j
          else:
            atom_begin = self.Z[j]
            atom_end = self.Z[i]
            index_begin = j
            index_end = i
          self.bonds[itr] = {'atom_begin'  : atom_begin,
                             'index_begin' : index_begin,
                             'atom_end'    : atom_end,
                             'index_end'   : index_end,
                             'length'      : d_ij}
          type_begin = ut.Z2n(atom_begin)
          type_end   = ut.Z2n(atom_end)
          bond_type  = type_begin + "-" + type_end
          if bond_type in self.bond_types:
            self.bond_types[bond_type] += 1
          else:
            self.bond_types[bond_type] = 1
          itr += 1

  def remove_atom(self, index):
    index -= 1
    if index < self.N - 1:
      out = copy.deepcopy(self)
      out.N -= 1
      out.R = np.delete(out.R, index, 0)
      out.Z = np.delete(out.Z, index)
      out.type_list = list(np.delete(out.type_list, index))
      return out
    else:
      print "index:%d out of range, nothing has happend" % index+1


  def have_bond(self, type_a, type_b):
    result = False
    if '0' not in self.bonds:
      self.find_bonds()
    if ut.n2Z(type_a) > ut.n2Z(type_b):
      atom_begin = ut.n2Z(type_b)
      atom_end = ut.n2Z(type_a)
    else:
      atom_begin = ut.n2Z(type_a)
      atom_end = ut.n2Z(type_b)
    for key in self.bonds:
      if self.bonds[key]['atom_begin'] == atom_begin and \
         self.bonds[key]['atom_end'] == atom_end:
        print self.bonds[key]['atom_begin'],
        print self.bonds[key]['atom_end']
        result = True
    return result
#  def have_bond(self, type_a, type_b):
#
#    result = False
#
#    na1 = ut.n2Z(type_a)
#    nb1 = ut.n2Z(type_b)
#    na2 = ut.n2Z(type_b)
#    nb2 = ut.n2Z(type_a)
#    def _qt2ob(qtmol):
#      mol = ob.OBMol()
#      new_atom = ob.OBAtom()
#      for atom in xrange(qtmol.N):
#        new_atom = mol.NewAtom()
#        new_atom.SetAtomicNum(qtmol.Z[atom])
#        new_atom.SetVector(qtmol.R[atom][0],
#                           qtmol.R[atom][1],
#                           qtmol.R[atom][2])
#      mol.ConnectTheDots()
#      #print "yo"
#      return mol
#      del new_atom
#      del mol
#
#    bond = ob.OBBond()
#    atom_a = ob.OBAtom()
#    atom_b = ob.OBAtom()
#    obmol = _qt2ob(self)
#
#    for i in range(obmol.NumBonds()):
#      bond = obmol.GetBond(i)
#      atom_a = bond.GetBeginAtom()
#      atom_b = bond.GetBeginAtom()
#      za = atom_a.GetAtomicNum()
#      zb = atom_b.GetAtomicNum()
#      if (za == na1 and zb == nb1) or (za == na2 and zb == nb2):
#        result = True
#        #print "(%d,%d) or (%d,%d)" % (na1,nb1,na2,nb2)
#    return result
#    del bond, atom_a, atom_b

  def center(self, center_coord):
    center_matrix = np.kron(
      np.transpose(np.atleast_2d(np.ones(self.N))),
      center_coord
    )
    self.R = self.R - center_coord

  def setCenterFrame(self, center_coord, frame_vector):
    print "not implemented yet"

  def setMultiplicity(self, m, **kargs):
    self.multiplicity = m
    if not ('forced' in kargs and kargs['forced']):
      if not (m % 2 != (np.sum(self.Z) + self.charge) % 2):
        ve = np.vectorize(ut.n2ve)
        nve = sum(ve(self.type_list)) - self.charge
        sys.exit("ERROR from geometry.py->" + \
                 "Molecule.setMultiplicity: " + \
                 "charge %d " % self.charge + \
                 "and multiplicity %d " % m + \
                 "with %d valence electrons " % nve +\
                 "are not compatible")

  def rotate(self, u, angle):
    print "not yet implemented"	

  def align(self, i,j,k):
    print "not yet implemented"

  def stretch(self):
    print "not yet implemented"

  def twist(self):
    print "not yet implemented"

  def sort(self):
    data = np.hstack([np.transpose(np.atleast_2d(self.Z)), self.R])
    data = data[data[:,0].argsort()]
    Z2n_vec = np.vectorize(ut.Z2n)
    self.Z = data[:,0].astype(int)
    self.type_list = Z2n_vec(self.Z)
    self.R = data[:,1:4]

    index_a = np.insert(self.Z, 0, 0)
    index_b = np.insert(self.Z, len(self.Z), 0)
    self.index = np.where((index_a != index_b))[0]

  def sort_coord(self, **kwargs):
    if 'order' in kwargs:
      order = kwargs['order']
    else:
      order = [0,1,2]
    ind = np.lexsort((self.R[:,order[2]],\
                      self.R[:,order[1]],\
                      self.R[:,order[0]]))
    self.R = self.R[ind]
    
    
#    self.Z = data[:,1:4]

  # read structrue from xyz
  def read_xyz(self, name, **kwargs):

    # caution! not format check. 
    # correct xyz format is assumed

    # local array varaible for easy append function
    coord = []
    type_list = []
    Z = []

    # open xyz file
    xyz_in = open(name, 'r')
    self.N = int(xyz_in.readline())
    xyz_in.readline()

    # loop through every line in xyz file
    for i in xrange(0, self.N):
      data = re.sub("[\n\t]", "",xyz_in.readline()).split(' ')
      # remove empty elements
      data = filter(None, data)
      type_list.append(data[0])
      Z.append(ut.n2Z(data[0]))
      crd = [float(data[1]),float(data[2]),float(data[3])]
      coord.append(crd)
    self.R = np.vstack(coord)
    self.type_list = np.array(type_list)
    self.Z = np.array(Z)

    if 'set_charge' in kwargs and kwargs['set_charge']:
      if np.sum(self.Z) % 2 == 1 :
        self.charge = -1


    xyz_in.close()

  # write xyz format to file
  def write_xyz(self, name):

    out=sys.stdout if re.match("stdout",name) else open(name,"w")

    #if len(self.type_list) != self.N:
    tlist = np.vectorize(ut.Z2n)
    self.type_list = tlist(self.Z)

    print >>out, str(self.N)+"\n"
    for I in xrange(0, self.N):
      print >>out, "%-2s " % self.type_list[I],
      print >>out, " ".join("% 8.4f" % i for i in self.R[I][:])

    if not re.match("stdout",name):
      out.close()

  # read structure from CPMD input
  def read_cpmdinp(self, name):
 
    self.N = 0
    #self.NType = 0
    #NTypeName = []
    coord = []
    Z = []
    type_list = []

    element_p = re.compile('\*([A-Za-z])*_')
    pp_p = re.compile('^\*')
    inp = open(name, 'r')
    while True:
      line = inp.readline()
      if not line: break
      if(re.match("&ATOMS",line)):
        while not (re.match("&END",line)):
          line = inp.readline()
          if(re.match(pp_p,line)):
            #self.NType += 1
            element = element_p.match(line).group(1)
            #NTypeName.append(element)
            inp.readline()
            N = int(inp.readline())
            for i in xrange(0, N):
              self.N += 1
              line = inp.readline()
              coord.append([float(x) for x in line.split()])
              type_list.append(element)
              Z.append(ut.n2Z(element))
    self.R = np.vstack(coord)
    #self.NTypeName = np.array(NTypeName)
    self.type_list = np.array(type_list)
    self.Z = np.array(Z)

    inp.close()

