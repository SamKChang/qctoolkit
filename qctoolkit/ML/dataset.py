import qctoolkit as qtk
import descriptors as dcr
from compiler.ast import flatten
import glob, re, os, copy
import random as rd
import numpy as np

# general object for a single data point
class DataPoint(object):
  def __init__(self, descriptor, *args, **kwargs):
    self._ref = np.nan
    self._input, self._name = getattr(dcr, descriptor)\
                                (*args, **kwargs).getVector()
  def setRef(ref):
    self._ref = ref

class DataSet(object):
  def __init__(self, path, pattern, descriptor, **kwargs):
    self.path = re.sub('/$', '', path)
    self.pattern = pattern
    self.data = []
    if 'threads' in kwargs:
      self.threads = int(kwargs['threads'])
    else:
      self.threads = 1

    if not os.path.exists(self.path):
      qtk.exit(self.path + "not found")
    qtk.report("DataSet", "reading folder", self.path,\
               "using descriptor:", descriptor)
    if descriptor == 'CoulombMatrix':
      if 'matrix_size' not in kwargs:
        qtk.warning("matrix size not assigend, " + \
                    "using default value")
        qtk.warning("matrix size WILL CHANGE " + \
                    "according to numer of atoms in the molecule")
        self.matrix_size = 0
      else:
        self.matrix_size = kwargs['matrix_size']

    if self.threads > 1:
      data_list = []
      for data in sorted(\
                    glob.glob(self.path + '/' + self.pattern)):
        data_list.append([descriptor, 
                          self.matrix_size, 
                          {'xyz':data}])
      self.data = qtk.parallelize(DataPoint,
                                  data_list,
                                  self.threads)
    else:
      for data in sorted(\
                    glob.glob(self.path + '/' + self.pattern)):
        self.data.append(\
          DataPoint(descriptor, self.matrix_size, xyz=data)\
        )

    self.data_size = len(self.data)

  def setReference(self, ref_vec):
    if len(ref_vec) != self.data_size:
      qtk.exit("number of data points not match")
    def set_ref(i, ref):
      self.data[i]._ref = ref
    vset_ref = np.vectorize(set_ref)
    vset_ref(range(self.data_size), ref_vec)
   

  def getTrainingSet(self, size):
    self.traning_size = size
    template = copy.deepcopy(self.data)
    index = range(self.data_size)
    rd.shuffle(index)
    index = index[:size]
    self.training_set = [template[i] for i in index]
    self.test_set = [data for i, 
                     data in enumerate(template) if i not in index]
  def prediction():
    pass
