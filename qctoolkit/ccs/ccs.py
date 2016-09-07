#from qctoolkit import *
import qctoolkit as qtk
import numpy as np
#from geometry import *
#from utilities import *
import re, copy, sys, os
from compiler.ast import flatten
import xml.etree.ElementTree as ET
import random
import yaml

class CCS(object):
  def __init__(self, xyz_file, parameter_file, **kwargs):
    self.structure = qtk.toMolecule(xyz_file)
    # mutation related variables
    self.mutation_list = []
    self.mutation_target = []
    self.mutation_size = 0
    # stretching related variables
    self.stretching_list = []
    self.stretching_direction = []
    self.stretching_range = []
    # rotation related variables
    self.rotation_list = []
    self.rotation_center = []
    self.rotation_axis = []
    self.rotation_range = []
    # replacing realted variables
    # not yet implemented

    # constraint variables #
    self.constraint = False
    self.forbiden_bonds = []
    self.ztotal = 0
    self.vtotal = 0
    self.element_count = {}

    # setup all parameters
    self.read_param(parameter_file)
    self.coor = flatten([
                 ['m' for _ in flatten(self.mutation_list)],
                 ['s' for _ in flatten(self.stretching_list)],
                 ['r' for _ in flatten(self.rotation_list)],
                ])

    MList = self.mutation_list
    _flatten = [item for sublist in MList for item in sublist]
    vlen = np.vectorize(len)
    try:
      lenList = vlen(MList)
    except TypeError:
      lenList = [len(MList[0]) for i in range(len(MList))]

    if not qtk.setting.quiet:
      report_itr = False
      if self.mutation_list:
        report_itr += True
        qtk.report('', "===== CCS REPORT =====", color=None)
        qtk.report("generating molecule", xyz_file)
        qtk.report("ccs parameter file", parameter_file)
        qtk.report("mutation indices", self.mutation_list)
        qtk.report("target atomic numbers", self.mutation_target)
        qtk.report("length of mutation vector",
               len(_flatten), "<=>", lenList)
        #print ""
      if self.stretching_list:
        qtk.report("stretching indices", self.stretching_list)
        qtk.report("stretching range", self.stretching_range)
        qtk.report("stretching direction indices",
               self.stretching_direction)
        #print ""
      if self.rotation_list:
        qtk.report("rotation indices", self.rotation_list)
        qtk.report("rotation center", self.rotation_center)
        qtk.report("rotation axis", self.rotation_axis)
        qtk.report("rotation range", self.rotation_range)
        #print ""
      qtk.status("ccs coordinate", self.coor)
      qtk.report('', "========= END ========", color=None)


  # !!!!! TODO !!!!! #
  # 100 line of read_param can be replace by simple xml reader
  # easier to maintain and extend

  def read_param(self,parameter_file):
    try:
      stem, extension = os.path.splitext(parameter_file)
      if extension == '.txt':
        self.read_param_txt(parameter_file)
      elif extension == '.xml':
        self.read_param_xml(parameter_file)
      elif extension == '.yml' or extension == '.yaml':
        self.read_param_yml(parameter_file)
      else:
        qtk.exit("extension " + extension + " not reconized...")
    except AttributeError:
      self.read_param_yml(parameter_file)
    size_list = [len(self.mutation_target[i]) ** \
                 len(self.mutation_list[i])\
                 for i in range(len(self.mutation_list))
                ]
    self.mutation_size = sum(size_list)

  ###########################################
  # convert data string to numerical values #
  ###########################################
  # used for xml/yaml format
  def str2data(self, data_string, **kwargs):
    if 'dtype' in kwargs:
      dtype = kwargs['dtype']
    else:
      dtype = 'index'
    _data = re.sub(re.compile('[ \n]*'),'',data_string)\
            .split(',')
    _out_list = []
    if dtype == 'index':
      for ind in _data:
        if re.match(re.compile(".*:.*"), ind):
          if re.match(re.compile("^[^0-9]*:"), ind):
            ind = re.sub("^[^0-9]*:", "0:", ind)
          if re.match(re.compile(":[^0-9]"), ind):
            ind = re.sub(":[^0-9]*", ":-1", ind)
          ind = map(int, ind.split(":"))
          for i in range(ind[0], ind[1]+1):
            _out_list.append(i)
        else:
          _out_list.append(int(ind))
    elif dtype == 'range':
      _out_list = map(float, _data[0].split(":"))
    return _out_list
  ##### end of data/string conversion #####

  #################################
  # read ccs_param from yaml file #
  #################################
  def read_param_yml(self, parameter_file):
    if type(parameter_file) is str:
      f = open(parameter_file, "r")
      param = yaml.safe_load(f)
    elif type(parameter_file) is dict:
      param = copy.deepcopy(parameter_file)

    # !!!!!!!!!!!!!!!!!
    # read span section
    def read_span(span):
      if span['type'] == 'mutation':
        self.mutation_list.append(self.str2data(span['index']))
        self.mutation_target.append(self.str2data(span['range']))
      elif span['type'] == 'stretching':
        self.stretching_list.append(self.str2data(span['index']))
        self.stretching_direction.append(self.str2data(\
          span['direction_index']))
        self.stretching_range.append(self.str2data(\
          span['range'], dtype='range'))
      else:
        qtk.exit(span['type'] + ' is not yet implemented')

    def read_constraint(constraint):
      if constraint['type'] == 'element_count':
        self.element_count.update(\
          {constraint['element']:\
           self.str2data(constraint['count'])})

    if param.has_key('span'):
      for span in param['span'].itervalues():
        read_span(span)
    if param.has_key('replace'):
      for replace in param['replace'].itervalues():
        pass # not yet implemented
    if param.has_key('constraint'):
      for constraint in param['constraint'].itervalues():
        read_constraint(constraint)

    if type(parameter_file) is str:
      f.close()
  ##### END OF READING YAML FILE #####

  ###############################################
  # read ccs_param to xml element tree directly #
  ###############################################
  def read_param_xml(self, parameter_file):
    tree = ET.parse(parameter_file)
    etroot = tree.getroot()

    self.mutation_size = 0
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # read span section of xml file
    def read_span(span):
      for _list in span:
        if _list.attrib['type']=='mutation':
          self.mutation_list.append(self.str2data(_list.text))
          self.mutation_target.append(self.str2data(\
            _list.attrib['range']))
        elif _list.attrib['type']=='stretching':
          self.stretching_list.append(self.str2data(_list.text))
          self.stretching_range.append(self.str2data(\
            _list.attrib['range'], dtype='range'))
          self.stretching_direction.append(self.str2data(\
            _list.attrib['fromto_direction_index']))
    # end of reading span section

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # read constraint section of xml file
    def read_constraint(constraint):
      for cst in constraint:
        if cst.tag == 'forbiden_bond':
          self.constraint = True
          self.forbiden_bonds.append(sorted(map(int, 
            [cst.attrib['begin'], cst.attrib['end']])))
        elif cst.tag == 'Z_sum':
          self.constraint = True
          self.ztotal = cst.attrib['total']
        elif cst.tag == 've_sum':
          self.constraint = True
          self.vtotal = cst.attrib['total']
        elif cst.tag == 'element_count':
          element_type = cst.attrib['type']
          element_count = int(cst.attrib['count'])
          self.element_count.update({element_type:element_count})
        else:
          qtk.exit("constraint '" + cst.tag + \
                   "' is not reconized")
    # end of reading constraint section

    for param in etroot:
      if param.tag == 'span':
        read_span(param)
      elif param.tag == 'constraint':
        self.constraint = True
        read_constraint(param)
  ##### END OF READING PARAM FROM XML #####

  ##########################################
  # GENERATE STRUCTURE FROM CCS COORDINATE #
  ##########################################
  # !!TODO!!
  # interface between mutation and other operations is necessary!
  def generate(self, **kwargs):
    self.new_structure = copy.deepcopy(self.structure)
    if 'suffix' in kwargs:
      self.new_structure.name = self.new_structure.name + "_" + \
                                kwargs['suffix']
      del kwargs['suffix']
    if 'mutation' in kwargs:
      self._mutate(kwargs['mutation'])
    if 'stretching' in kwargs:
      self._stretch(kwargs['stretching'])
    if 'rotation' in kwargs:
      self._rotate(kwargs['rotatting'])
    return self.new_structure

  def _mutate(self, mutation):
    for m in xrange(len(mutation)):
      for i in xrange(len(mutation[m])):
        index = self.mutation_list[m][i] - 1
        #target = self.mutation_target[m][mutation[m][i]]
        target = mutation[m][i]
        self.new_structure.Z[index] = target
        self.new_structure.type_list[index] = qtk.Z2n(target)
  def _stretch(self, stretching):
    pass
  def _rotate(self, rotation):
    pass
  ##### END OF STRUCTURE GENERATION #####

  ################################################
  # GENERATOR TO CONSTRUCT NEXT STRUCTURE IN CCS #
  ################################################
  def iterateMutation(self, **kwargs):
    if 'space' not in kwargs:
      space = 1
    else:
      space = kwargs['space']
    digits = len(str(self.mutation_size / space))
    itr = 0
    i = 0
    while i < self.mutation_size:
      mutation_crd = self.mutationNumber(i)
      crd = self.random()[1]
      crd['mutation'] = mutation_crd
      crd['suffix'] = str(itr).rjust(digits, '0')
      i = i + space
      itr = itr + 1
      yield self.generate(**crd)
  ##### END OF NEXT STRUCTURE #####

  #################################################
  # INDEX CONVERSION FROM INTEGER TO CCS MUTATION #
  #################################################
  def mutationNumber(self, query_int):
    out = []
    constant = 0
    too_large = True
    for i in range(len(self.mutation_list)):
      ilist = self.mutation_list[i]
      itarg = self.mutation_target[i]
      base = len(itarg)
      imax = base ** len(ilist)
      icrd = [0 for j in range(len(ilist))]
      if query_int < imax:
        number = query_int - constant
        ind = qtk.numberToBase(number, base)
        icrd[:len(ind)] = ind
        too_large = False
      else:
        constant = constant + imax
        icrd = [base - 1 for j in range(len(ilist))]
      mcrd = [self.mutation_target[i][j] for j in icrd]
      out.append(mcrd)
    if not too_large:
      return out
    else:
      return None
  ##### END OF INDEX CONVERSION #####

  ################################################
  # TEST A GIVEN STRUCTURE SATISFIED CONSTRAINTS #
  ################################################
  # NOT YET IMPLEMENTED!!!!!
  def onManifold(self,structure):
    on_manifold = True
    if self.constraint:
#      if self.ztotal:
#        if self.ztotal == sum(structure.Z):
#          on_manifold = on_manifold*True
#        else:
#          on_manifold = on_manifold*False
#      if self.vtotal:
#        if self.vtotal == structure.getValenceElectrons():
#          on_manifold = on_manifold*True
#        else:
#          on_manifold = on_manifold*False
      return on_manifold
    else:
      return True

  ###########################################
  # GENERATE A RANDOM POINT ON CCS MANIFOLD #
  ###########################################
  def random(self, **kwargs):

    def map2list(index, nested_list):
      new_list = flatten(nested_list)
      value = new_list[index]
      for i in range(len(nested_list)):
        if value in nested_list[i]:
          j = nested_list[i].index(value)
          return i, j
      

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # generate a random ccs coordinate sequence
    def random_coord(self, mode):
      # determin coordinate type
      if mode == 'mutation':
        _dtype = int
        _list = copy.deepcopy(self.mutation_list)
        _targ = copy.deepcopy(self.mutation_target)


      # algorithm for element_count constraint
      # 1. construct mapping from flattened to nested mutation list
      # 2. select all candidate indices for constrainted element
      # 3. remove constrained element from target lists
      # 4. random choice of candidate for constrained mutation
      # 5. loop through all other indices without constrained elem

        # single constraint implementation
        if self.element_count:
          new_list = []
          new_targ = []
          constraint = []
          selected = []
          for elem in self.element_count.iterkeys():
            itr = 0
            count = 0
            _candidate = []
            for sublist in _targ:
              if qtk.n2Z(elem) in sublist:
                _candidate.extend(range(itr, 
                                        itr+len(_list[count])))
                _targ[count].remove(qtk.n2Z(elem))
              itr += len(_list[count])
              count += 1
            _candidate = [x for x in _candidate\
                          if not x in flatten(selected)]
            
            selection = sorted(random.sample(_candidate, 
                        random.choice(self.element_count[elem])),
                        reverse=True)
            selected.append(selection)
  
            for ind in selection:
              i, j = map2list(ind, _list)
              constraint.append([(i,j), qtk.n2Z(elem)])
          constraint_ind = [a[0] for a in constraint]
        # end of constraint setting

      elif mode == 'stretching':
        _list = self.stretching_list
        _targ = self.stretching_range
        _dtype = float
      elif mode == 'rotation':
        _list = self.rotation_list
        _targ = self.rotation_range
        _dtype = float

      # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      # main routine to generate a random coordinate
      _vec = []
      ind = 0
      for _group, _targp in zip(_list, _targ):
        _subvec = []
        for _atom in _group:
          if _dtype == int:
            # element_count constraint
            if self.element_count:
              i, j = map2list(ind, _list)
              if (i,j) in constraint_ind:
                index = constraint_ind.index((i,j))
                _subvec.append(constraint[index][1])
              else:
                _subvec.append(random.choice(_targp))
            # without constraint
            else:
              _subvec.append(random.choice(_targp))
          elif _dtype == float:
            _subvec.append(random.uniform(_targp[0], _targp[1]))
          ind += 1
        _vec.append(_subvec)
      if len(_vec) > 0:
        return _vec
    # end of random coordinate generation

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!
    # loop until a point is found
    while True:
      mut_vec = None
      while mut_vec is None:
        try:
          mut_vec = random_coord(self,'mutation')
        except:
          pass
      str_vec = random_coord(self,'stretching')
      rot_vec = random_coord(self,'rotation')
      query = {}
      if mut_vec:
        query.update({'mutation':mut_vec})
      if str_vec:
        query.update({'stretching':str_vec})

      new_structure = self.generate(**query)
      if self.onManifold(new_structure):
        break
    # end of loop, i.e. a ccs point is found

    return new_structure, query
  ##### END OF RANDOM CCS POINT #####

#  def addr2loc(self, inp, mode='mutation'):
#    if mode == 'mutation':
#      _list = self.mutation_list
#    if type(inp) is int:
#      elm = 0
#      grp = 0
#      for i in _list:
#        if inp < elm + len(i):
#          return grp, inp-elm
#        else:
#          elm = elm + len(i)
#          grp = grp + 1
#    if type(inp) is list:
      

  #########################################
  # mating function for genetic algorithm #
  #########################################
  def mate(self, parent1, parent2, mutation_rate):

    def getChild():
      child = {}
      for key in parent1.iterkeys():
        if key == 'mutation':
          # implementation for element_count constraints
          if self.element_count:
            constraint_list = []
            constraint_keys = []
            # tool to constuct list
            def _list_append(mcoord, i_list, Zc):
              for g in range(len(mcoord)):
                grp = mcoord[g]
                for i in range(len(grp)):
                  Z_p = grp[i]
                  if Z_p == Zc:
                    i_list.append((g,i))
            # construct mutation list
            for elem in self.element_count.iterkeys():
              i_list = []
              _Z = qtk.n2Z(elem)
              _list_append(parent1['mutation'], i_list, _Z)
              _list_append(parent2['mutation'], i_list, _Z)
              constraint_list.append(i_list)
              constraint_keys.append(elem)

            while True:
              selected_list = []
              selected_coord = {}
              consistant = True
              for i in range(len(constraint_list)):
                ind_location = constraint_list[i]
                random.shuffle(ind_location)
                elem = qtk.n2Z(constraint_keys[i])
                _count = self.element_count[constraint_keys[i]][0]
                _end = _count
                while len(set(ind_location[:_end])) < _count:
                  _end += 1
                _new = list(set(ind_location[:_end]))
                to_select = copy.deepcopy(selected_list)
                to_select.extend(_new)
                if len(set(to_select)) == \
                len(selected_list) + len(_new):
                  selected_list.extend(_new)
                  for coord in _new:
                    selected_coord.update({coord:elem})
                  consistant = consistant & True
                else:
                  consistant = False
              child_ind = []
              for g in range(len(parent1['mutation'])):
                child_list = []
                for i in range(len(parent1['mutation'][g])):
                  if selected_coord.has_key((g,i)):
                    child_list.append(selected_coord[(g,i)])
                  else:
                    candidate = []
                    A1 = qtk.Z2n(parent1['mutation'][g][i])
                    A2 = qtk.Z2n(parent2['mutation'][g][i])
                    if not self.element_count.has_key(A1):
                      candidate.append(qtk.n2Z(A1))
                    if not self.element_count.has_key(A2):
                      candidate.append(qtk.n2Z(A2))
                    if len(candidate) == 0:
                      consistant = False
                    else:
                      child_list.append(random.choice(candidate))
                child_ind.append(child_list)
              if consistant:
                break
            for g in range(len(parent1['mutation'])):
              for i in range(len(parent1['mutation'][g])):
                if not selected_coord.has_key((g,i)):
                  if random.random() < mutation_rate:
                    child_ind[g][i] = random.choice(\
                      self.mutation_target[g])
                    elem = qtk.Z2n(child_ind[g][i])
                    if self.element_count.has_key(elem):
                      coord = random.choice([tupl for tupl,targ 
                        in selected_coord.iteritems()\
                        if targ == qtk.n2Z(elem)])
                      del selected_coord[coord]
                      selected_coord.update({(g,i):qtk.n2Z(elem)})
                      [gc, ic] = list(coord)
                      new_targ = random.choice([z for z in \
                                  self.mutation_target[gc]\
                                  if not self.element_count\
                                  .has_key(qtk.Z2n(z))])
                      child_ind[gc][ic] = new_targ
            child.update({'mutation':child_ind})
              
          else:
            child_mut = []
            for i in range(len(parent1[key])):
              # mix two parent
              grp = parent1[key][i]
              index = range(len(grp))
              random.shuffle(index)
              ind1 = index[:len(index)/2]
              if len(index)%2 == 0:
                ind2 = index[len(index)/2:]
                append = -1
              else:
                ind2 = index[len(index)/2:-1]
                append = index[-1]
              new = [parent1[key][i][ind] for ind in ind1]
              new.extend([parent2[key][i][ind] for ind in ind2])
              if append >= 0:
                if random.random() >= 0.5:
                  new.append(parent1[key][i][append])
                else:
                  new.append(parent2[key][i][append])
              # mutation, loop through each element
              for j in range(len(new)):
                if random.random() < mutation_rate:
                  mut_target = [tar \
                    for tar in self.mutation_target[i]\
                    if tar != new[j]]
                  new[j] = random.choice(mut_target)
              child_mut.append(new)
            child.update({'mutation':child_mut})
        else:
          qtk.exit("ccs mate function is not yet implemented for "\
                   + key)
      return child

    child = getChild()
    child_structure = self.generate(**child)
    while not self.onManifold(child_structure):
      child = getChild()
      child_structure = self.generate(**child)
    return child
  

  #########################################
  # read ccs_param from plain text format #
  #########################################
  # deprecated. kept for backward sporting
  def read_param_txt(self, parameter_file):
    param = open(parameter_file, 'r')

    mutation_flag = re.compile(".*mutation_list:.*")
    stretching_flag = re.compile(".*stretching_list:.*")
    direction_flag = re.compile(".*direction_index:.*")
    rotation_flag = re.compile(".*rotation_list:.*")
    center_flag = re.compile(".*center_index:.*")
    end_flag = re.compile(".*end.*")

    while True:
      line = param.readline()
      if not line: break
      line = re.sub(re.compile("#.*"),'', line)
      if re.match(mutation_flag, line.lower()):
        while not re.match(end_flag, line):
          line = param.readline().rstrip()
          line = re.sub(re.compile("#.*"),'', line)
          try:
            mlist_str =  re.sub(re.compile("->.*"), "", line)\
                         .split(",")
            tlist_str =  re.sub(re.compile(".*->"), "", line)\
                         .split(",")
            mlist = []
            tlist = []
            for m in mlist_str:
              if re.match(re.compile(".*:.*"), m):
                rangeList = map(int, m.split(":"))
                mmin = rangeList[0]
                mmax = rangeList[1] + 1
                mlist.extend(range(mmin,mmax))
              else:
                mlist.append(int(m))
            for t in tlist_str:
              if re.match(re.compile(".*:.*"), t):
                rangeList = map(int, t.split(":"))
                tmin = rangeList[0]
                tmax = rangeList[1] + 1
                tlist.extend(range(tmin,tmax))
              else:
                tlist.append(int(t))
            self.mutation_list.append(mlist)
            self.mutation_target.append(tlist)
          except ValueError:
            pass
        MList = self.mutation_list
        _flatten = [item for sublist in MList for item in sublist]
        vlen = np.vectorize(len)
        # wired numpy.int32 bug for equal length sublist
        try:
          lenList = vlen(MList)
        except TypeError:
          lenList = [len(MList[0]) for i in range(len(MList))]
      elif re.match(stretching_flag, line.lower()):
        while not re.match(end_flag, line):
          line = param.readline().rstrip()
          line = re.sub(re.compile("#.*"),'', line)
          if re.match(direction_flag, line):
            dlist = map(int, re.sub(re.compile(".*:"),
                                    "", line).split(","))
            self.stretching_direction.append(dlist)
            line = param.readline().rstrip()
            line = re.sub(re.compile("#.*"),'', line)
            try:
              slist = map(int, re.sub(re.compile("->.*"),\
                                      "", line).split(","))
              self.stretching_list.append(slist)
              rlist = map(float, re.sub(re.compile(".*->"),\
                                        "", line).split(':'))
              self.stretching_range.append(rlist)
            except ValueError:
              pass

      elif re.match(rotation_flag, line.lower()):
        while not re.match(end_flag, line):
          line = param.readline()
          line = re.sub(re.compile("#.*"),'', line)
          if re.match(center_flag, line):
            center = int(re.sub(re.compile(".*:"),
                                "", line))
            self.rotation_center.append(center)
            line = param.readline().rstrip()
            line = re.sub(re.compile("#.*"),'', line)
            axis = map(int, re.sub(re.compile(".*:"),
                                   "", line).split(","))
            self.rotation_axis.append(axis)
            line = param.readline().rstrip()
            line = re.sub(re.compile("#.*"),'', line)
            try:
              rolist = map(int, re.sub(re.compile("->.*"),\
                                       "", line).split(","))
              self.rotation_list.append(rolist)
              rglist = map(float, re.sub(re.compile(".*->"),\
                                         "", line).split(':'))
              self.rotation_range.append(rglist)
            except ValueError:
              pass
    param.close()

