from geometry import *
import re, copy

class MoleculeSpan(object):
  def __init__(self, xyz_file, parameter_file):
    self.structure = Molecule()
    self.structure.read_xyz(xyz_file)
    # mutation related variables
    self.mutation_list = []
    self.mutation_target = []
    # stretching related variables
    self.stretching_direction = []
    self.stretching_list = []
    self.stretching_range = []
    # rotation related variables
    self.rotation_center = []
    self.rotation_axis = []
    self.rotation_list = []
    self.rotation_range = []
    self.read_param(parameter_file)

  def read_param(self, parameter_file):
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
      if re.match(mutation_flag, line.lower()):
        while not re.match(end_flag, line):
          line = param.readline().rstrip()
          try:
            mlist = map(int, re.sub(re.compile("->.*"), "", line)\
                            .split(","))
            tlist = map(int, re.sub(re.compile(".*->"), "", line)\
                             .split(","))
            self.mutation_list.append(mlist)
            self.mutation_target.append(tlist)
          except ValueError:
            pass
      elif re.match(stretching_flag, line.lower()):
        while not re.match(end_flag, line):
          line = param.readline().rstrip()
          if re.match(direction_flag, line):
            dlist = map(int, re.sub(re.compile(".*:"),
                                    "", line).split(","))
            self.stretching_direction.append(dlist)
            line = param.readline().rstrip()
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
          if re.match(center_flag, line):
            center = int(re.sub(re.compile(".*:"),
                                "", line))
            self.rotation_center.append(center)
            line = param.readline().rstrip()
            axis = map(int, re.sub(re.compile(".*:"),
                                   "", line).split(","))
            self.rotation_axis.append(axis)
            line = param.readline().rstrip()
            try:
              rolist = map(int, re.sub(re.compile("->.*"),\
                                       "", line).split(","))
              self.rotation_list.append(rolist)
              rglist = map(float, re.sub(re.compile(".*->"),\
                                         "", line).split(':'))
              self.rotation_range.append(rglist)
            except ValueError:
              pass

  def generate(self, **kwargs):
    if 'mutation' in kwargs:
      return self._mutate(kwargs['mutation'])
    if 'stretching' in kwargs:
      self._stretch(kwargs['stretching'])
    if 'rotation' in kwargs:
      self._rotate(kwargs['rotatting'])

  def _mutate(self, mutation):
    out = copy.deepcopy(self.structure)
    for m in xrange(len(mutation)):
      for i in xrange(len(mutation[m])):
        index = self.mutation_list[m][i] - 1
        target = self.mutation_target[m][mutation[m][i]]
        out.Z[index] = target
    return out
  def _stretch(self, stretching):
    print stretching
  def _rotate(self, rotation):
    print rotation
