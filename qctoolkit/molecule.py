import numpy as np
import qctoolkit as qtk
import setting
import re, os, sys, copy, operator
from time import sleep
import networkx
from networkx.algorithms.components.connected\
     import connected_components
import periodictable as pt
import collections
from math import ceil, floor

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
      alignSVD(mol, ref_list=[], tar_list=[]) --- align to mol according
                                                  to SVD minimization
                                                  default ref_list = 
                                                          tar_list = all
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
  def __init__(self, *args, **kwargs):
    # number of atoms
    self.N = 0
    # atom coordinates
    self.R = np.atleast_2d(np.array([]))
    self.R_scale = np.atleast_2d(np.array([]))
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
    self.symmetry = False
    self.grid = False
    self.name = ''

    if len(args) == 1:
      mol = args[0]
      self.read(mol, **kwargs)
    elif len(args) == 2:
      N = len(args[0])
      dim1 = np.array(args[0]).shape
      dim2 = np.array(args[1]).shape
      if dim1 == (N,) and dim2 == (N, 3):
        atoms = args[0]
        coord = np.array(args[1])
      elif dim1 == (N, 3) and dim2 == (N,):
        atoms = args[1]
        coord = np.array(args[0])
      else:
        qtk.exit('not supported declaration of molecule object.')
      self.addAtoms(atoms, coord)
     

    attr_list = dir(self)
    for string, value in kwargs.iteritems():
      if string in attr_list:
        setattr(self, string, kwargs[string])

    if self.N > 0:
      self.ve = self.getValenceElectrons()
      self.ne = sum(self.Z)

  def __repr__(self):
    if self.name: return self.name
    else: return 'generic Molecule object'

  # tested
  def __add__(self, other):
    if self.R_scale.shape[1] > 0:
      qtk.exit('Molecule add not implemented for crystals.' + \
               'use extend/setAtoms instead')
    out = Molecule()
    out.N = self.N + other.N
    out.R = np.vstack([self.R, other.R])
    out.Z = np.hstack([self.Z, other.Z])
    out.type_list = np.hstack([self.type_list, other.type_list])
    out.string = np.hstack([self.string, other.string])
    out.charge = self.charge + other.charge
    out.name = self.name + "_" + other.name
    return out

  def copy(self):
    return copy.deepcopy(self)

  def nuclear_repulsion(self):
    out = 0.
    for i in range(self.N):
      for j in range(i+1, self.N):
        Rij = np.linalg.norm(self.R[i] - self.R[j]) * 1.8897261245650618
        out -= self.Z[i] * self.Z[j] / Rij
    return out

  def view(self, name=None):
    tmp = copy.deepcopy(self)
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
    if 'output' not in kwargs:
      kwargs['output'] = 'string'
    if kwargs['output'] == 'dictionary' or kwargs['output'] == 'count':
      out = {}
      for element in data:
        out[qtk.Z2n(element[0])] = element[1]
    elif kwargs['output'] == 'string':
      out = ''
      for element in data:
        out = out + qtk.Z2n(element[0]) + str(element[1])
    return out

  def build(self, moleculeData=None, name=None, **kwargs):
    if moleculeData is not None:
      if type(moleculeData) is list:
        moleculeData = np.array(moleculeData)
      if len(moleculeData.shape) == 1:
        moleculeData = np.atleast_2d(moleculeData)
      self.N = moleculeData.shape[0]
      self.Z = moleculeData[:, 0]
      self.R = moleculeData[:, 1:]
      self.type_list = [qtk.Z2n(z) for z in self.Z]
      self.string = ['' for i in range(self.N)]
      if name is None:
        self.name = self.stoichiometry()
      else:
        self.name = name
    else:
      assert 'Z' in kwargs
      assert 'R' in kwargs
      Z = np.array(kwargs['Z'])
      R = np.array(kwargs['R'])
      self.N = len(Z)
      self.Z = Z
      self.R = R
      self.type_list = [qtk.Z2n(z) for z in self.Z]
      self.string = ['' for i in range(self.N)]
      self.name = self.stoichiometry()
    return self

  # tested
  def findBonds(self, ratio=setting.bond_ratio, **kwargs):
    del self.segments
    del self.bond_types
    self.segments = []
    self.bond_types = {}
    if 'no_report' not in kwargs or not kwargs['no_report']:
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
          bond_table = qtk.data.elements.bond_table
          bond_keys = []
          bond_keys = [
            type_begin + _ + type_end for _ in ['-', '=', '#']
          ]
          try:
            bond_type_ind = np.argmin(
              abs(
                np.array([
                          bond_table[k][0] for k in bond_keys
                          if k in bond_table.keys()
                         ]) - d_ij
              )
            )
          except Exception as _e:
            self.write_xyz()
            qtk.exit(
              "error while processing bond" +\
              str(bond_keys) + "with error message %s" % str(_e))
          bond_type = bond_keys[bond_type_ind]
          self.bonds[itr]['name'] = bond_type
          bond_energy = \
            bond_table[bond_keys[bond_type_ind]][1] * \
            qtk.convE(1, 'kj-kcal')[0]
          self.bonds[itr]['energy'] = bond_energy
          if np.isnan(bond_energy):
            qtk.warning("Non-tabliated covalent bond %s" % bond_type)
          if bond_type in self.bond_types:
            self.bond_types[bond_type] += 1
          else:
            self.bond_types[bond_type] = 1
          itr += 1

    segments = list(connected_components(to_graph(bond_list)))
    for s in range(len(segments)):
      segment = list(segments[s])
      new_mol = self.getSegment(segment, **kwargs)
      ns = len(self.segments)
      new_mol.name = new_mol.name + '_%d' % ns
      self.segments.append(new_mol)
    for s in [i for i in range(self.N) if not bonded[i]]:
      segment = [s]
      new_mol = self.getSegment(segment, **kwargs)
      ns = len(self.segments)
      new_mol.name = new_mol.name + '_%d' % ns
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
    type_a = type_a.title()
    type_b = type_b.title()
    if '0' not in self.bonds:
      self.findBonds()
    if qtk.n2Z(type_a) > qtk.n2Z(type_b):
      atom_begin = qtk.n2Z(type_b)
      atom_end = qtk.n2Z(type_a)
    else:
      atom_begin = qtk.n2Z(type_a)
      atom_end = qtk.n2Z(type_b)
    for key in self.bonds:
      if self.bonds[key]['atom_begin'] == atom_begin \
      and self.bonds[key]['atom_end'] == atom_end:
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
  def center(self, center_coord = None):
    if center_coord is None:
      center_coord = self.getCenterOfMass()
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
    shift = np.kron(vector, template)
    self.R += shift

  # tested
  def align(self, u=None, **kwargs):
    if 'no_center' in kwargs and kwargs['no_center']:
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
      if np.linalg.norm(u[1]) - u[1] > 1E-8 or u[0] > 0:
        u = np.array(u)
        n1 = np.linalg.norm(u)
        n2 = np.linalg.norm(v)
        angle = np.arccos(np.dot(u, v)/(n1*n2))
      else:
        u = np.array([0, 0, 1])
        angle = np.pi

      #if angle > 0:
      axis = np.cross(u, v)
      axis = axis / np.linalg.norm(axis)
      self.rotate(angle, axis)

  def alignSVD(self, mol, ref_list=None, tar_list=None):
    if type(mol) is str:
      try:
        mol = qtk.Molecule(mol)
      except:
        qtk.exit("error when reading molecule file: %s" % mol)
    assert issubclass(mol.__class__, qtk.Molecule)
    if not ref_list:
      ref_list = [i for i in range(self.N)]
    if not tar_list:
      tar_list = copy.deepcopy(ref_list)

    lst_a = self.R[ref_list]
    lst_b = mol.R[tar_list]
    center_a = np.mean(lst_a, axis=0)
    center_b = np.mean(lst_b, axis=0)
    #na = len(lst_a)
    na = self.N
    #nb = len(lst_b)
    nb = mol.N
    crd_a = self.R - np.kron(center_a, np.ones((self.N, 1)))
    crd_b = mol.R - np.kron(center_b, np.ones((mol.N, 1)))
    ref_a = lst_a - np.kron(center_a, np.ones((len(lst_a), 1)))
    ref_b = lst_b - np.kron(center_a, np.ones((len(lst_b), 1)))

    H = np.dot(np.transpose(ref_a), ref_b)
    U, s, V = np.linalg.svd(H)
    R = np.dot(np.transpose(V), np.transpose(U))
    self.R = np.transpose(
      np.dot(R, np.transpose(crd_a))) + \
      np.kron(center_b, np.ones((na, 1))
    )

  def alignAtoms(self, ind1, ind2, ind3):
    self.center(self.R[ind1])
    vec = self.R[ind2] - self.R[ind1]
    self.align(vec)
    self.center(self.R[ind1])

    if abs(self.R[ind3][1]) > 0:
      tangent = self.R[ind3][2]/self.R[ind3][1]
    else:
      if self.R[ind3][2] > 0:
        tangent = np.inf
      else:
        tangent = -np.inf

    angle = np.arctan(tangent)
    self.rotate(-angle, [1,0,0])
    self.center(self.R[ind1])
    if self.R[ind3][1] < 0:
      self.rotate(np.pi, [1, 0, 0])
      

  def rotate(self, angle, u):
    R_tmp = copy.deepcopy(self.R)
    self.R = np.dot(qtk.R(angle, u),R_tmp.T).T

  def twist(self, targets, vec_ind, angle_degree):
    angle = angle_degree / 180.0 * np.pi
    center = copy.deepcopy(self.R[vec_ind[0]])
    self.center(self.R[vec_ind[0]])
    R_part = copy.deepcopy(self.R[targets])
    vector = self.R[vec_ind[1]] - self.R[vec_ind[0]]
    vector = vector / np.linalg.norm(vector)
    R_part = np.dot(qtk.R(angle, vector), R_part.T).T
    self.R[targets] = R_part
    self.shift(center)

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
    if type(element) is not list:
      element = [element]
      coord = [coord]
    if type(element[0]) is not str:
      element = [qtk.Z2n(int(i)) for i in element]
    Z = list(self.Z)
    for i in range(len(element)):
      #e = getattr(pt, element[i].title())
      e = qtk.element[element[i]]
      r = coord[i]
      if self.N == 0:
        self.R = np.array(r)
      else:
        self.R = np.vstack([self.R, np.array(r)])
      self.N = self.N + 1
      self.type_list.append(e.symbol)
      Z.append(e.number)
      self.string.append('')
    self.Z = np.array(Z)
    if self.N == 1:
      self.R = np.array([self.R])

  # tested
  def setAtoms(self, index, **kwargs):
    if type(index) is int:
      index = [index]
    if 'element' in kwargs or 'Z' in kwargs:
      for i in index:
        if 'element' in kwargs:
          if type(kwargs['element']) is str:
            Z = qtk.n2Z(kwargs['element'])
            Zn = kwargs['element']
          elif type(kwargs['element']) is int\
          or type(kwargs['element']) is float:
            Z = kwargs['element']
            Zn = qtk.Z2n(Z)
        elif 'Z' in kwargs:
          Z = kwargs['Z']
          Zn = qtk.Z2n(Z)
        self.Z[i] = Z
        self.type_list[i] = Zn
    if 'string' in kwargs:
      minZ = min(min(self.Z)-1, 0)
      for i in index:
        self.string[i] = kwargs['string']
        self.Z[i] = minZ

  def mutateElement(self, oldZ, newZ):
    if oldZ != newZ:
      if type(oldZ) is str:
        oldZ = qtk.n2Z(oldZ)
      if type(newZ) is str:
        newZ = qtk.n2Z(newZ)
      newElem = qtk.Z2n(newZ)
      indices = []
      for i in range(self.N):
        if abs(self.Z[i] - oldZ) < 0.1:
          indices.append(i)
      self.setAtoms(indices, element=newElem)

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

    self.grid = [[mol_min[i], mol_max[i], step[i]] for i in range(3)]
    return self.grid

  def setCelldm(self, celldm=None, **kwargs):
    if not self.celldm:
      if not celldm:
        if not self.periodic:
          if 'margin' not in kwargs:
            margin = qtk.setting.box_margin
          else:
            margin = kwargs['margin']
          self.box = self.getSize() + \
                     2*np.array([margin, margin, margin])
          self.celldm = copy.deepcopy(self.box)
          self.celldm.extend([0, 0, 0])
      else:
        self.periodic = True
        self.celldm = celldm
        self.box = celldm[:3]
        self.R_scale = qtk.xyz2fractional(self.R, celldm)
    else:
      if celldm:
        self.periodic = True
        for i in range(len(celldm)):
          self.celldm[i] = celldm[i]
        self.box = celldm[:3]
        for i in range(self.N):
          for j in range(3):
            if self.scale:
              self.R[i, j] = self.R_scale[i, j] * \
                             self.celldm[j] / float(self.scale[j])
            else:
              self.R[i, j] = self.R_scale[i, j] * celldm[j]
    return self.celldm

  def expand(self, ratio):
    if self.celldm:
      for i in range(3):
        self.celldm[i] = self.celldm[i] * ratio
      self.R = qtk.fractional2xyz(self.R_scale, self.celldm)
    else:
      self.R = self.R * ratio

  def extend(self, ratio):

    def take(data, mask):
      return list(data[i] for i in range(len(mask)) if mask[i])

    assert len(ratio) == 3
    assert len(self.R_scale) == self.N
    factor = reduce(lambda x, y: x*y, ratio)
    max_R = [ceil(i) for i in np.max(self.R_scale, axis = 0)]
    for i in range(3):
      M = self.N
      R_scale = copy.deepcopy(self.R_scale)
      R = copy.deepcopy(self.R)
      Z = np.array(self.Z)
      type_list = np.array(self.type_list)
      string = np.array(self.string)
      mask = [True for j in range(M)]
      for r in range(int(floor(ratio[i])) - 1):
        S = M
        new_R_scale = copy.deepcopy(R_scale)
        new_R = copy.deepcopy(R)
        for j in range(len(new_R_scale)):
          new_scale_j = new_R_scale[:, i][j] + max_R[i] * (r + 1)
          new_Rj = new_R[:, i][j] + self.celldm[i] * (r + 1)
          if new_scale_j < max_R[i] * ratio[i]:
            new_R_scale[:, i][j] = new_scale_j
            new_R[:, i][j] = new_Rj
          else:
            S = S - 1 
            mask[j] = False
        self.N = S + self.N
        self.R_scale = np.vstack(
                         [self.R_scale, 
                          take(new_R_scale,mask)]
                       )
        self.R = np.vstack([self.R, take(new_R, mask)])
        self.Z = list(np.hstack([self.Z, take(Z,mask)]))
        new_list = np.hstack([self.type_list, take(type_list,mask)])
        self.type_list = [str(a) for a in new_list]
        new_str = np.hstack([self.string, take(string,mask)])
        self.string = [str(s) for s in new_str]
      self.celldm[i] = self.celldm[i] * ratio[i]
      self.scale =  [ceil(i) for i in np.max(self.R_scale, axis = 0)]

  def copy(self):
    return copy.deepcopy(self)

  # tested by qminp
  def sort(self, order = 'Zxyz'):
    odict = {'x':0, 'y':1, 'z':2}
    tmp = []
    for o in order:
      if o == 'Z':
        tmp.insert(0, self.Z)
      elif o in odict:
        tmp.insert(0, self.R[:, odict[o]])
      else:
        qtk.exit("sorting order '%c' not valid" % o)
    ind = np.lexsort(tmp)
    self.R = self.R[ind]
    if list(self.R_scale[0]):
      self.R_scale = self.R_scale[ind]
    self.Z = list(np.array(self.Z)[ind])
    self.type_list = list(np.array(self.type_list)[ind])
    self.string = list(np.array(self.string)[ind])

    if order == 'Zxyz':
      type_list = []
      self.index = []
      for i in range(self.N):
        Zn = self.type_list[i]
        if Zn not in type_list:
          type_list.append(Zn)
          self.index.append(i)
      self.index.append(self.N)
    
  def sort_coord(self, **kwargs):
    if 'order' in kwargs:
      order = kwargs['order']
    else:
      order = 'xyz'
    self.sort(order)

  # tested
  # general interface to dertermine file type
  def read(self, name, **kwargs):
    if os.path.exists(name):
      stem, extension = os.path.splitext(name)
      if re.match('\.xyz', extension):
        self.read_xyz(name, **kwargs)
      elif re.match('\.ascii', extension):
        self.read_ascii(name, **kwargs)
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

    prop_list = ['charge', 'celldm', 'scale', 'symmetry']
    for prop in prop_list:
      try:
        prop_str = filter(lambda x: prop in x, content)[0]
        prop_str = re.sub('.*:', '', prop_str)
        prop_data = prop_str.split(' ')
        prop_data = filter(None, prop_data)
        if len(prop_data) == 1:
          try:
            setattr(self, prop, float(prop_data[0]))
          except Exception as exc:
            if prop == 'symmetry':
              setattr(self, prop, prop_data[0].strip())
        elif len(prop_data) > 1:
          setattr(self, prop, [float(_) for _ in prop_data])
      except ValueError as exc:
        setattr(self, prop, False)
        qtk.warning("setting attribute %s with error: %s" % \
          (prop, exc))
      except:
        setattr(self, prop, False)
    
    if not self.charge:
      self.charge = 0
    if self.celldm or self.scale:
      self.periodic = True
    if self.celldm:
       assert len(self.celldm) == 6

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

      self.R_scale = copy.deepcopy(self.R)
      self.R = qtk.fractional2xyz(self.R_scale, self.celldm)

  def read_ascii(self, name, **kwargs):
    def getParameters(content):
      line1 = content[1]
      line2 = content[2]
      # v_sim cell format
      # a_x b_x b_y
      # c_x c_y c_z
      [a, b_cosC, b_sinC] = [
        float(p) for p in filter(None, line1.split(' '))
      ]
      [c_cosB, c_cosAsinC, c_z] = [
        float(p) for p in filter(None, line2.split(' '))
      ]
      b = np.sqrt(b_cosC**2 + b_sinC**2)
      gamma = np.arccos(b_cosC / b)
      c = np.sqrt(c_z**2 + c_cosB**2 + c_cosAsinC**2)
      beta = np.arccos(c_cosB / c)
      alpha = np.arccos(c_cosAsinC / c / np.sin(gamma))
      return [a, b, c, np.cos(alpha), np.cos(beta), np.cos(gamma)]

    def getMolecule(content):
      comment_p = re.compile(r'^ *\t*(?!([!#])).*$')
      data_p = re.compile(r'^[\ \t0-9\-\. ]+.*$')
      dataContent = filter(comment_p.match, content[3:])
      dataContent = filter(data_p.match, dataContent)
      N = len(dataContent)
      R = []
      molData = []
      for line in dataContent:
        data = filter(None, line.replace('\n', '').split(' '))
        molLine = [qtk.n2Z(data[3])]
        molLine.extend([float(r) for r in data[:3]])
        molData.append(molLine)
      return molData

    xyz = open(name, 'r')
    content = xyz.readlines()
    xyz.close()
    self.build(getMolecule(content))
    celldm = getParameters(content)
    self.setCelldm(celldm)

  # tested
  def write(self, *args, **kwargs):
    if len(args) > 0:
      _, ext = os.path.splitext(args[0])
      kwargs['format'] = ext[1:]
    if 'format' not in kwargs:
      kwargs['format'] = 'xyz'
    if kwargs['format'] == 'xyz':
      self.write_xyz(*args, **kwargs)
    elif kwargs['format'] == 'pdb':
      self.write_pdb(*args, **kwargs)
    elif kwargs['format'] == 'ascii':
      self.write_ascii(*args, **kwargs)
    #elif kwargs['format'] == 'cyl':
    #  self.write_cyl(*args, **kwargs)
    else:
      qtk.exit("output format: " + kwargs['format'] \
               + " not reconized")

  # tested
  # write xyz format to file
  def write_xyz(self, name=None, **kwargs):

    def listStr(data_list):
      out = str(data_list)
      out = out.replace('[', '')
      out = out.replace(']', '')
      return out.replace(',', '')

    out = sys.stdout if not name else open(name,"w")

    if 'fractional' in kwargs and kwargs['fractional']:
      if self.celldm and self.scale:
        out.write(str(self.N)+"\n\n")
        for I in xrange(0, self.N):
          out.write("%-2s " % self.type_list[I])
          out.write(" ".join("% 8.4f" % i \
            for i in self.R_scale[I][:]))
          out.write('\n')
        out.write("\ncelldm: %s\n" % listStr(self.celldm))
        out.write('scale: %s\n' % listStr(self.scale))
        if self.symmetry:
          out.write('symmetry: %s\n' % self.symmetry)
      else:
        del kwargs['fractional']
        qtk.warning('celldm or scale not set, print cartician')
        self.write(name, **kwargs)
    else:
      out.write(str(self.N)+"\n\n")
      for I in xrange(0, self.N):
        out.write("%-2s " % self.type_list[I])
        out.write(" ".join("% 8.4f" % i for i in self.R[I][:]))
        out.write("\n")
  
      if name:
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
    out.write("%-10s%s\n" % ('COMPND', self.stoichiometry()))
    out.write("%-10s%s\n" % ('AUTHOR', 'QCTOOLKIT'))
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
        out.write("%-6s%5d%3s%6s%6d     " %\
          ('ATOM', itr+i, atom.upper(), 'LIG', chain) +\
          "% 7.3f % 7.3f % 7.3f%6.2f%6.2f%12s\n" %\
          (xi, yi, zi, 1, 0, atom))
      connection = connect(segment, itr, connection)
      itr = itr + segment.N
      chain = chain + 1
    print >> out, connection,
    print >> out, "END"
    if not re.match("",name):
      out.close()

  def write_ascii(self, name = None, **kwargs):
    assert len(self.celldm) == 6
    cdm = self.celldm
    gamma = np.arccos(cdm[5])
    ax = cdm[0]
    bx = cdm[1] * cdm[5]
    by = cdm[1] * np.sin(gamma)
    cx = cdm[2] * cdm[4]
    cy = cdm[2] * cdm[3] * np.sin(gamma)
    cz = np.sqrt(cdm[2]**2 - cy**2 - cx**2)

    out = sys.stdout if not name else open(name,"w")
    out.write("%d\n" % self.N)
    out.write("% 9.6f % 9.6f % 9.6f\n" % (ax, bx, by))
    out.write("% 9.6f % 9.6f % 9.6f\n" % (cx, cy, cz))
    for I in xrange(0, self.N):
      out.write(" ".join("% 8.4f" % i for i in self.R[I][:]))
      out.write(" %-2s " % self.type_list[I])
      out.write("\n")
    out.close()
