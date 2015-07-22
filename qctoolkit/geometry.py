import numpy as np
import utilities as ut
import re
import sys

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

  def __add__(self, other):
    out = Molecule()
    out.N = self.N + other.N
    out.R = np.vstack([self.R, other.R])
    out.Z = np.hstack([self.Z, other.Z])
    return out

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
    
    
#    self.Z = data[:,1:4]

  # read structrue from xyz
  def read_xyz(self, name):

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

