import numpy as np
#import utilities as ut
import qctoolkit as qtk
import setting
import re, os, sys, copy, operator
from time import sleep
import networkx
from networkx.algorithms.components.connected import connected_components
import periodictable as pt
import collections

class Molecule(object):
  """
  Molecule class for basic molecule operation

  args: 
    mol(str), molecule file

  attributes: 
    N(int), number of atoms
    R(np.array), coordinate array
    Z(list(int/float)), list of nuclear charges
    type_list
    charge
    multiplicity
    bonds
    bond_types
    segments
    name

  methods:
    connectivity/property:
      stoichiometry(format=)
      findBonds(ratio=None) --- get all connected atom pairs within
                                cutoff distance = atom radius*ratio
                                where ratio is set to default 1.1
      haveBond('A', 'B') --- return True/False depends on the 
                             existance of connection A-B
                             where A, B are element symbols
      setChargeMultiplicity(c, m) --- set charge multiplicity to c, m
      setCharge(c) --- set charge to c, under same multiplicity
      getValenceElectrons() --- return number of valence electrons
    geometry operation:
      distance(i, j) --- return distance (angstrom) between 
                         atom i and j
      center(coord) --- shift molecule such that coord becomes 
                        zero vector. It can be used for centering 
                        selected atom: R = R - coord
      shift(vector) --- shift molecule by vector: R = R + vector
      align(vector, axis) --- align vecotr to axis, axis is set to
                              [1,0,0] by default. It takes 0=x, 1=y,
                              2=z, or a 3D vector
      stretch([i,j], [s,t], d) --- stretch atom i,j along direction of
                                   atom s,t to distance d
      rotate() ---
      twist() ---
      getCenter() --- return 3x1 array of center of coordinates
      getCenterOfCharge() --- return 3x1 array of center of charge
      getCenterOfMass() --- return 3x1 array of center of mass
      principalAxes() --- return eigenvalue/vectors of 
                          momentum of inertia tensor
    modify molecule:
      addAtoms(list_str, list_coord) --- add atoms in list_str with 
                                         coordinates list_coord
      removeAtoms(list_ind) --- remove list of atoms in list_ind
      setAtoms(list_ind, list_str) --- set list of atoms to strings
                                       in list_str. Used for user 
                                       define atom symbols
      isolateAtoms(list_ind) --- keep only the atoms listed in list_ind
    basic IO:
      read()
      read_xyz()
      read_pdb() --- not yet implemented
      write()
      write_xyz()
      write_pdb()

    QM related:
      getSize() --- return 1x3 array where each component denotes
                    the space spanned by the molecule
  """

  # used for pymol numeration
  mol_id = 0
  def __init__(self, mol=None, **kwargs):
    # number of atoms
    self.N = 0
    # atom coordinates
    self.R = np.atleast_2d(np.array([]))
    # atom symbols
    self.type_list = []
    # nuclear charges
    self.Z = []
    # moelcule charge
    self.charge = 0
    self.multiplicity = 1
    # index of different atoms
    self.index = 0
    self.bonds = {}
    self.bond_types = {}
    self.string = []
    self.segments = []
    self.periodic = False
    self.scale = False
    self.celldm = False
    self.grid = False
    self.name = ''
    if mol:
      self.read(mol, **kwargs)

  def __repr__(self):
    if self.name: return self.name
    else: return 'generic Molecule object'

  # tested
  def __add__(self, other):
    out = Molecule()
    out.N = self.N + other.N
    out.R = np.vstack([self.R, other.R])
    out.Z = np.hstack([self.Z, other.Z])
    out.type_list = np.hstack([self.type_list, other.type_list])
    out.string = np.hstack([self.string, other.string])
    out.charge = self.charge + other.charge
    out.name = self.name + "_" + other.name
    return out

  def view(self, name=None):
    tmp = copy.deepcopy(self)
    if self.scale:
      try:
        for i in range(3):
          tmp.R[:,i] = tmp.R[:,i] * tmp.celldm[i]\
                       / float(tmp.scale[i])
      except AttributeError:
        pass

    if qtk.imported('pymol'):
      qtk.report("Molecule", "initializing pymol...", color=None)
      import pymol
      pymol.finish_launching()
    else:
      pymol.cmd.reinitialize()
      sleep(0.5)
    if name:
      tmp_file = name + "_tmp_" + str(Molecule.mol_id) + '.xyz'
    else:
      tmp_file = 'pymol_tmp_' + str(Molecule.mol_id) + '.xyz'
    Molecule.mol_id = Molecule.mol_id + 1
    tmp.write_xyz(tmp_file)
    pymol.cmd.load(tmp_file)
    os.remove(tmp_file)

  # tested
  def stoichiometry(self, **kwargs):
    elements = collections.Counter(sorted(self.Z))
    data = zip(elements.keys(), elements.values())
    data.sort(key=lambda tup: tup[0])
    if 'format' not in kwargs:
      kwargs['format'] = 'string'
    if kwargs['format'] == 'list':
      return data
    elif kwargs['format'] == 'string':
      out = ''
      for element in data:
        out = out + qtk.Z2n(element[0]) + str(element[1])
      return out

  # tested
  def findBonds(self, ratio=setting.bond_ratio, **kwargs):
    del self.segments
    self.segments = []
    qtk.report("Molecule", 
               "finding bonds with cutoff ratio", 
               ratio)
    def to_graph(l):
      G = networkx.Graph()
      for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(to_edges(part))
      return G
    
    def to_edges(l):
      """ 
      treat `l` as a Graph and returns it's edges 
      to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
      """
      it = iter(l)
      last = next(it)
    
      for current in it:
        yield last, current
        last = current 
    itr = 0
    bond_list = []
    bonded = [False for i in range(self.N)]
    for i in xrange(self.N):
      for j in xrange(i+1, self.N):
        d_ij = np.linalg.norm(self.R[i,:] - self.R[j,:])
        atom_i = getattr(pt, self.type_list[i])
        atom_j = getattr(pt, self.type_list[j])
        Ri = atom_i.covalent_radius + \
             atom_i.covalent_radius_uncertainty
        Rj = atom_j.covalent_radius + \
             atom_j.covalent_radius_uncertainty
        Dij = (Ri+Rj) * float(ratio)
        if d_ij < Dij:
          bonded[i] = True
          bonded[j] = True
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
          bond_list.append([i, j])
          type_begin = qtk.Z2n(atom_begin)
          type_end   = qtk.Z2n(atom_end)
          bond_type  = type_begin + "-" + type_end
          if bond_type in self.bond_types:
            self.bond_types[bond_type] += 1
          else:
            self.bond_types[bond_type] = 1
          itr += 1

    segments = list(connected_components(to_graph(bond_list)))
    for s in range(len(segments)):
      segment = list(segments[s])
      new_mol = self.getSegment(segment, **kwargs)
      self.segments.append(new_mol)
    for s in [i for i in range(self.N) if not bonded[i]]:
      segment = [s]
      new_mol = self.getSegment(segment, **kwargs)
      self.segments.append(new_mol)

  # tested
  def getSegment(self, index_list, **kwargs):
    new_mol = copy.deepcopy(self)
    new_mol.charge = 0
    if type(index_list) != list:
      index_list = [index_list]
    index_list = map(lambda a: a, index_list)
    new_mol.N = len(index_list)
    new_mol.R = np.array([new_mol.R[i] for i in index_list])
    new_mol.Z = np.array([new_mol.Z[i] for i in index_list])
    new_mol.type_list = \
      [new_mol.type_list[i] for i in index_list]
    new_mol.string = np.array(new_mol.string)[index_list].tolist()
    unpaired = new_mol.getValenceElectrons() % 2
    if unpaired == 1:
      if 'charge_saturation' not in kwargs:
        new_mol.setChargeMultiplicity(-1, 1, **kwargs)
      else:
        assert type(kwargs['charge_saturation']) is int
        charge = kwargs['charge_saturation']
        new_mol.setChargeMultiplicity(charge, 1, **kwargs)
    return new_mol

  # tested
  def haveBond(self, type_a, type_b):
    result = False
    if '0' not in self.bonds:
      self.findBonds()
    if qtk.n2Z(type_a) > qtk.n2Z(type_b):
      atom_begin = qtk.n2Z(type_b)
      atom_end = qtk.n2Z(type_a)
    else:
      atom_begin = qtk.n2Z(type_a)
      atom_end = qtk.n2Z(type_b)
    for key in self.bonds:
      if self.bonds[key]['atom_begin'] == atom_begin and \
         self.bonds[key]['atom_end'] == atom_end:
        print self.bonds[key]['atom_begin'],
        print self.bonds[key]['atom_end']
        result = True
    return result

  # tested
  def getValenceElectrons(self):
    ve = np.vectorize(qtk.n2ve)
    nve = sum(ve(self.type_list)) - self.charge
    return nve

  # tested
  def setChargeMultiplicity(self, c, m, **kwargs):
    if type(c) == int or type(c) == float:
      self.charge = c
    if type(m) == int:
      self.multiplicity = m

    if type(self.multiplicity)==int and\
       type(self.charge)==(int or float):
      if not (self.multiplicity % 2 != \
              (np.sum(self.Z) + self.charge) % 2):
        ve = np.vectorize(qtk.n2ve)
        nve = sum(ve(self.type_list)) - self.charge
        msg = "Multiplicity %d " % self.multiplicity + \
              "and %d valence electrons " % nve +\
              "\n(with charge %3.1f) " % float(self.charge) +\
              "are not compatible"
        qtk.prompt(msg + "\nsuppress warning py no_warning=True,"\
                  + " continue?")

  # tested
  def setCharge(self, **kwargs):
    self.charge = 0
    unpaired = self.getValenceElectrons() % 2
    if unpaired == 1:
      if 'charge_saturation' not in kwargs:
        self.setChargeMultiplicity(-1, 1)
      else:
        assert type(kwargs['charge_saturation']) is int
        charge = kwargs['charge_saturation']
        self.setChargeMultiplicity(charge, 1)

  # tested
  def distance(self, i, j):
    Ri = self.R[i]
    Rj = self.R[j]
    return np.linalg.norm(Ri - Rj)

  # tested
  def center(self, center_coord):
    center_matrix = np.kron(
      np.transpose(np.atleast_2d(np.ones(self.N))),
      center_coord
    )
    self.R = self.R - center_coord

  # tested
  def shift(self, shift_vector):
    shift_matrix = np.kron(
      np.transpose(np.atleast_2d(np.ones(self.N))),
                   np.array(shift_vector)
    )
    self.R = self.R + shift_matrix

  # tested
  def stretch(self, targets, direction_indices, distance):
    if type(targets) is not list:
      targets = [targets]
    direction = [self.R[index] for index in direction_indices]
    vector = direction[1] - direction[0]
    vector = distance * vector/np.linalg.norm(vector) + direction[0]
    template = np.zeros([self.N,1])
    template[targets] = 1
    ref = copy.deepcopy(self.R)
    ref[targets,:] = 0
    ref = self.R - ref
    shift = np.kron(vector, template) - ref
    self.R += shift

  # tested
  def align(self, u=None, **kwargs):
    center = self.getCenterOfMass()
    self.center(center)
    if 'axis' not in kwargs:
      v = np.array([1,0,0])
    else:
      if type(kwargs['axis']) is int:
        assert kwargs['axis'] >= 0 and kwargs['axis'] < 3
        if kwargs['axis'] == 0:
          v = np.array([1,0,0])
        elif kwargs['axis'] == 1:
          v = np.array([0,1,0])
        elif kwargs['axis'] == 2:
          v = np.array([0,0,1])
      else:
        v = np.array(kwargs['axis'])
        assert v.shape == (3,)

    if u is None:
      U = self.principalAxes()[1]
      self.R = np.dot(self.R, U)
    else:      
      u = np.array(u)
      n1 = np.linalg.norm(u)
      n2 = np.linalg.norm(v)
      angle = np.arccos(np.dot(u, v)/(n1*n2))
      axis = np.cross(u, v)
      axis = axis / np.linalg.norm(axis)
      self.rotate(angle, axis)

  def rotate(self, angle, u):
    R_tmp = copy.deepcopy(self.R)
    self.R = np.dot(qtk.R(angle, u),R_tmp.T).T

  def twist(self):
    print "not yet implemented"

  # tested
  def getCenter(self):
    return np.sum(self.R, axis=0)/self.N

  # tested
  def getCenterOfCharge(self):
    weighted = self.R * np.array(self.Z).reshape([self.N,1])
    return np.sum(weighted, axis=0)/float(sum(self.Z))

  # tested
  def getCenterOfMass(self):
    mass_list = [qtk.n2m(elem) for elem in self.type_list]
    weighted = self.R * np.array(mass_list).reshape([self.N,1])
    return np.sum(weighted, axis=0)/float(sum(mass_list))

  # tested
  def principalAxes(self, **kwargs):
    weight = [qtk.n2m(elem) for elem in self.type_list]
    center = self.getCenterOfMass()
    self.center(center)

    inertia = np.zeros([3,3])
    I0 = 0
    for i in range(3):
      I0 = I0 + sum(self.R[:,i]**2 * weight)
    for i in range(3):
      coord_i = self.R[:,i]
      inertia[i,i] = I0 - sum(coord_i**2 * weight)
      for j in range(i+1,3):
        coord_j = self.R[:,j]
        inertia[i,j] = -sum(coord_i*coord_j*weight)
        inertia[j,i] = inertia[i,j]
    I, U = np.linalg.eigh(inertia)
    self.shift(center)
    #return sorted(I,reverse=True), U[I.argsort()[::-1]]
    return I, U

  # tested
  def addAtoms(self, element, coord):
    def getAtom(element):
      if len(element)<2:
        try:
          atom = getattr(pt, element.title())
          return atom
        except:
          qtk.exit("element %s not found." % element.title())
      else:
        try:
          atom = getattr(pt, element.lower())
          return atom
        except:
          qtk.exit("element %s not found." % element.lower())
    
    if type(element) is str:
      element = [element]
      coord = [coord]
    Z = list(self.Z)
    for i in range(len(element)):
      e = getAtom(element[i])
      r = coord[i]
      if self.N == 0:
        self.R = np.array(r)
      else:
        self.R = np.vstack([self.R, np.array(r)])
      self.N = self.N + 1
      self.type_list.append(e.symbol)
      Z.append(e.number)
    self.Z = np.array(Z)

  # tested
  def setAtoms(self, index, **kwargs):
    newZ = self.Z
    if type(index) is int:
      index = [index]
    if 'element' in kwargs:
      for i in index:
        assert i>=0
        if 'element' in kwargs:
          Z = qtk.n2Z(kwargs['element'])
        elif 'Z' in kwargs:
          Z = kwargs['Z']
        newZ[i] = Z
        self.type_list[i] = qtk.Z2n(Z)
    if 'string' in kwargs:
      minZ = min(min(self.Z)-1, 0)
      for i in index:
        self.string[i] = kwargs['string']
        newZ[i] = minZ
    self.Z = newZ

  # tested
  def removeAtoms(self, indices):
    if type(indices) is int:
      indices = [indices]
    indices.sort()
    i_max = indices[-1]
    if i_max <= self.N - 1:
      for i in range(len(indices)):
        index = indices[len(indices) - 1 - i]
        self.N -= 1
        self.R = np.delete(self.R, index, 0)
        self.Z = np.delete(self.Z, index)
        self.type_list = list(np.delete(self.type_list, index))
    else:
      msg = "index:%d out of range:%d, nothing has happend"\
            % (index, self.N)
      qtk.warning(msg)

  # tested
  def isolateAtoms(self, index_list, **kwargs):
    if type(index_list) != list:
      index_list = [index_list]
    self.N = len(index_list)
    self.R = np.array([self.R[i] for i in index_list])
    self.Z = np.array([self.Z[i] for i in index_list])
    self.type_list = \
      [self.type_list[i] for i in index_list]
    self.string = np.array(self.string)[index_list].tolist()

  # tested by qminp
  def getSize(self):
    """
    return box dimension spanned by the coordinates
    """
    def size(i):
      return max(self.R[:,i]) - min(self.R[:,i])
    return np.array([size(i) for i in range(3)])

  def setGrid(self, margin=None, **kwargs):
    if not margin: 
      margin = qtk.setting.box_margin
    mol_max = [max(self.R[:,i]) + margin for i in range(3)]
    mol_min = [min(self.R[:,i]) - margin for i in range(3)]
    size = self.getSize() + 2*np.array([margin, margin, margin])

    if 'step' not in kwargs:
      step = [int(size[i] / 0.12) for i in range(3)]
    else:
      step = kwargs['step']

    print type(self)
    print self.grid
    self.grid = [[mol_min[i], mol_max[i], step[i]] for i in range(3)]
    print self.grid
    return self.grid

  def setCelldm(self, margin=None):
    if not self.periodic:
      if not margin:
        margin = qtk.setting.box_margin
      self.box = self.getSize() + 2*np.array([margin, margin, margin])
      self.celldm = copy.deepcopy(self.box)
      self.celldm.extend([0, 0, 0])
    return self.celldm

  # tested by qminp
  def sort(self):
    new = sorted(zip(self.R, self.type_list, self.Z, self.string), 
                 key=operator.itemgetter(2))
    self.R = np.array([_R for _R, _T, _Z, _S in new])
    self.type_list = np.array([_T for _R, _T, _Z, _S in new])
    self.Z = np.array([_Z for _R, _T, _Z, _S in new])
    self.string = np.array([_S for _R, _T, _Z, _S in new])
   
    index_a = np.insert(self.Z, 0, 0)
    index_b = np.insert(self.Z, len(self.Z), 0)
    self.index = np.where((index_a != index_b))[0]
    if self.index[0] != 0:
      self.index = np.insert(self.index, 0, 0)
    if index_a[-1] == index_b[-1] and index_a[-1] == 0:
      self.index = sorted(np.insert(self.index, \
                                    0, len(self.index)))
  # tested by qminp
  def sort_coord(self, **kwargs):
    if 'order' in kwargs:
      order = kwargs['order']
    else:
      order = [0,1,2]
    ind = np.lexsort((self.R[:,order[2]],\
                      self.R[:,order[1]],\
                      self.R[:,order[0]]))
    self.R = self.R[ind]

  # tested
  # general interface to dertermine file type
  def read(self, name, **kwargs):
    if os.path.exists(name):
      stem, extension = os.path.splitext(name)
      if re.match('\.xyz', extension):
        self.read_xyz(name, **kwargs)
      else:
        qtk.exit("suffix " + extension + " is not reconized")

      self.string = ['' for _ in range(self.N)]
      self.name = qtk.fileStrip(stem)

      if np.sum(self.Z) % 2 == 1:
        if 'charge_saturation' not in kwargs:
          self.charge = -1
        else: 
          self.charge = kwargs['charge_saturation']
    else:
      qtk.exit("file: '" + name + "' not found")

  def read_xyz(self, name, **kwargs):
    xyz = open(name, 'r')
    content = xyz.readlines()
    xyz.close()
    content = [line.replace('\t', ' ') for line in content]

    prop_list = ['charge', 'celldm', 'scale']
    for prop in prop_list:
      try:
        prop_str = filter(lambda x: prop in x, content)[0]
        prop_str = re.sub('.*:', '', prop_str)
        prop_data = prop_str.split(' ')
        if len(prop_data) == 1:
          setattr(self, prop, float(prop_data[0]))
        elif len(prop_data) > 1:
          setattr(self, prop, [float(_) for _ in prop_data])
      except:
        setattr(self, prop, False)
    if not self.charge: self.charge = 0
    if self.celldm or self.scale: self.periodic = True

    self.N = int(content[0])
    coord_list = content[2 : self.N + 2]
    coord = [filter(None,[a for a in entry.split(' ')]) 
             for entry in coord_list]
    type_list = list(np.array(coord)[:,0])
    self.type_list = [str(elem) for elem in type_list]
    self.Z = [qtk.n2Z(elem) for elem in self.type_list]
    self.Z = np.array(self.Z)
    self.R = np.array(coord)[:,1:4].astype(float)

    self.box = False
    if self.celldm:
      self.periodic = True
      angle = self.celldm[3:]
      angle_sum = sum([abs(entry) for entry in angle])
      if angle_sum == 0:
        self.box = self.celldm[:3]
        
      if self.scale:
        if angle_sum == 0:
          self.R_scale = self.R
          print np.array(self.scale)
          factor = np.array(self.box) / np.array(self.scale)
          factor = np.kron(factor, np.ones((self.N, 1)))
          self.R = self.R * factor
        else:
          msg = 'non cubic box with scale is not supported'
          qtk.warning(msg)

  # tested
  def write(self, *args, **kwargs):
    if 'format' not in kwargs:
      kwargs['format'] = 'xyz'
    if kwargs['format'] == 'xyz':
      self.write_xyz(*args, **kwargs)
    elif kwargs['format'] == 'pdb':
      self.write_pdb(*args, **kwargs)
    #elif kwargs['format'] == 'cyl':
    #  self.write_cyl(*args, **kwargs)
    else:
      qtk.exit("output format: " + kwargs['format'] \
               + " not reconized")

  # tested
  # write xyz format to file
  def write_xyz(self, *args, **kwargs):
    if len(args) == 1:
      name = args[0]
    else: name = ''
    out = sys.stdout if not name else open(name,"w")

    #if len(self.type_list) != self.N:
    #tlist = np.vectorize(qtk.Z2n)
    #self.type_list = list(tlist(self.Z))

    print >>out, str(self.N)+"\n"
    for I in xrange(0, self.N):
      print >>out, "%-2s " % self.type_list[I],
      print >>out, " ".join("% 8.4f" % i for i in self.R[I][:])

    if not re.match("",name):
      out.close()

  # tested
  # write pdb format to file
  # amino acid is not implemented!
  def write_pdb(self, *args, **kwargs):
    if len(args) == 1:
      name = args[0]
    else: 
      name = ''
    out = sys.stdout if not name else open(name,"w")
    if len(self.segments) == 0:
      self.findBonds(quiet=True)
    print >> out, "%-10s%s" % ('COMPND', self.stoichiometry())
    print >> out, "%-10s%s" % ('AUTHOR', 'QCTOOLKIT')
    chain = 1
    itr = 1

    def connect(molecule, shift, connection):
      molecule.findBonds(quiet=True)
      for i in molecule.bonds.iterkeys():
        bond = molecule.bonds[i]
        ai = bond['index_begin'] + shift
        aj = bond['index_end'] + shift
        connection = connection +\
          "%-6s%4d %4d\n" % ('CONECT', ai, aj)
      return connection

    connection = ''
    for segment in self.segments:
      for i in range(segment.N):
        atom = segment.type_list[i]
        xi = segment.R[i, 0]
        yi = segment.R[i, 1]
        zi = segment.R[i, 2]
        #"% 7.3f % 7.3f % 7.3f%6.2f%6.2f%12s" %\
        print >> out, "%-6s%5d%3s%6s%6d     " %\
          ('ATOM', itr+i, atom.upper(), 'LIG', chain) +\
          "% 7.3f % 7.3f % 7.3f%6.2f%6.2f%12s" %\
          (xi, yi, zi, 1, 0, atom)
      connection = connect(segment, itr, connection)
      itr = itr + segment.N
      chain = chain + 1
    print >> out, connection,
    print >> out, "END"
    if not re.match("",name):
      out.close()
