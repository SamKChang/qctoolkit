import qctoolkit as qtk
import descriptors as dcr
import kernels as ker
from compiler.ast import flatten
import glob, re, os, copy
import random as rd
import numpy as np
import kernel_matrix as km
import kernel_vectors as kv
import multiprocessing as mp
import operator

# general object for a single data point
class DataPoint(object):
  def __init__(self, descriptor, *args, **kwargs):
    self._vector, self._name = getattr(dcr, descriptor)\
                                (*args, **kwargs).getVector()
    self.ref = np.nan
    self.alpha = np.nan
    self.prediction = np.nan
  def getName(self):
    return self._name
  def getRef(self):
    return self.ref
  def getVector(self):
    return self._vector
  def getPrediction(self):
    return self.prediction

class DataSet(object):
  def __init__(self, path, pattern, **kwargs):
    self.path = re.sub('/$', '', path)
    self.pattern = pattern
    if not os.path.exists(self.path):
      qtk.exit(self.path + "not found")

  #################
  # som utilities #
  #################
  def getNameList(self):
    def get_name(i):
      return self.data[i].getName()
    vget_name = np.vectorize(get_name)
    return vget_name(range(self.data_size))

  def getRefList(self):
    def get_ref(i):
      return self.data[i].getRef()
    vget_ref = np.vectorize(get_ref)
    return vget_ref(range(self.data_size))

  ##################
  # SET DESCRIPTOR #
  ##################
  def setDescriptor(self, descriptor, **kwargs):
    self.descriptor = descriptor
    if 'threads' in kwargs:
      self.threads = int(kwargs['threads'])
    else:
      self.threads = 1
    qtk.report("DataSet", "reading folder", self.path)
    qtk.report("Descriptor", self.descriptor)
    if descriptor == 'CoulombMatrix':
      if 'matrix_size' not in kwargs:
        qtk.warning("matrix size not assigend, " + \
                    "using default value")
        qtk.warning("matrix size WILL CHANGE " + \
                    "according to numer of atoms in the molecule")
        self.matrix_size = 0
      else:
        self.matrix_size = kwargs['matrix_size']
    else:
      qtk.exit("descriptor" + descriptor + "is not implemented")

    if self.threads > 1:
      data_list = []
      for data in sorted(\
                    glob.glob(self.path + '/' + self.pattern)):
        data_list.append([descriptor, 
                          self.matrix_size, 
                          {'xyz':data}])
      self.data = qtk.parallelize(DataPoint,
                                  data_list,
                                  threads=self.threads)
    else:
      for data in sorted(\
                    glob.glob(self.path + '/' + self.pattern)):
        self.data.append(\
          DataPoint(descriptor, self.matrix_size, xyz=data)\
        )

    self.data_size = len(self.data)
  ###### END OF SET DESCRIPTOR ######

  ##############
  # SET KERNEL #
  ##############
  def setKernel(self, kernelName, *klargs):
    self.kernel = kernelName
    self.klargs = klargs
    qtk.report("Kernel", self.kernel ,list(self.klargs))
    #self.kernel = getattr(ker, kernelName)(*klargs, **klkwargs)

  ###### END OF SET KERNEL ######

  #################
  # SET REFERENCE #
  #################
  def setReference(self, ref_vec):
    if len(ref_vec) != self.data_size:
      qtk.exit("number of data points not match")
    def set_ref(i, ref):
      self.data[i].ref = ref
    vset_ref = np.vectorize(set_ref)
    vset_ref(range(self.data_size), ref_vec)
  ###### END OF SET REFERENCE ######
   
  ############
  # TRAINING #
  ############
  def training(self, size, **kwargs):
    """
    set up training set, test set, and calculate alphas
    """

    # 'lambda' is a python keyword...
    if 'deviation' in kwargs:
      self._lambda = kwargs['deviation']
    else:
      self._lambda = 0

    qtk.report("Deviation parameter", self._lambda)

    template = copy.deepcopy(self.data)
    index = range(self.data_size)
    rd.shuffle(index)
    max_index = size
    i = 0
    self.training_set = []
    self.training_index = []
    reference_vector = []
    ref_coord = []
    # keep the flexibility for unset reference datapoint
    while i<max_index:
      if i==self.data_size:
        qtk.exit("reference not set")
      data_point = template[index[i]]
      if not np.isnan(data_point.ref):
        self.training_set.append(data_point)
        self.training_index.append(index[i])
        reference_vector.append(data_point.getRef())
        ref_coord.append(data_point.getVector())
      else:
        max_index += 1
      i += 1
    self.refVector = np.array(reference_vector)
    self.refCoord = np.array(ref_coord)
    
    self.rest_set = [data for i, data in enumerate(template) 
                     if i not in self.training_index]
    def reset_alpha(i):
      self.rest_set[i].alpha = np.nan
    vreset_alpha = np.vectorize(reset_alpha)
    vreset_alpha(range(len(self.rest_set)))

    rows, columns = self.refCoord.shape

    qtk.progress("Kernel", "generating",\
                 "%dx%d" % (size, size), "kernel matrix...")
    # exteral C module
    self.kernelMatrix = km.kernel_matrix(self.refCoord,
                                         rows, columns,
                                         self.kernel,self.klargs)
    qtk.done()

    if 'eigen' in kwargs and kwargs['eigen']:
      qtk.progress("Kernel", "proceed diagonalization...")
      self.kernelEW, self.kernelEV =\
        np.linalg.eigh(self.kernelMatrix)
      qtk.done()

    qtk.progress("Kernel", "inverting...")
    self.alphas = np.dot(np.linalg.inv(self.kernelMatrix\
                         + self._lambda*np.eye(size))\
                        ,self.refVector)
    qtk.done()
    def set_alpha(i):
      self.training_set[i].alpha = self.alphas[i]
    vset_alpha = np.vectorize(set_alpha)
    vset_alpha(range(len(self.training_set)))
  ###### END OF TRAINING ######

  ###########
  # PREDICT #
  ###########
  def predict(self, size, **kwargs):
    template = copy.deepcopy(self.rest_set)
    index = range(len(self.rest_set))
    rd.shuffle(index)
    test_list = index[:size]
    i = 0
    self.test_set = [self.rest_set[i] for i in test_list]
    self.kernelVectors = np.atleast_2d([])
    self.testCoord = np.array(
                      [data.getVector() for data in self.test_set])
    rrows, rcolumns = self.refCoord.shape
    trows, tcolumns = self.testCoord.shape

    # exteral C module
    qtk.progress("Predition", "generating",\
                 "%dx%d" % (trows, rrows),\
                 "kernel projection..." )
    self.kernelVectors = kv.kernel_vectors(self.refCoord, 
                                           rrows, rcolumns, 
                                           self.testCoord, 
                                           trows, tcolumns,
                                           self.kernel,self.klargs)
    qtk.done()
    test_true = []
    test_pred = []
    for i in range(len(self.test_set)):
      prediction = np.dot(self.kernelVectors[i], self.alphas)
      self.test_set[i].prediction = prediction
      test_true.append(self.test_set[i].ref)
      test_pred.append(prediction)

    # analize and status report
    self.testTrue = np.array(test_true)
    self.testPred = np.array(test_pred)
    self.error = abs(self.testPred - self.testTrue)
    self.MAE = sum(self.error)/len(self.error)
    self.RMSE = np.sqrt(sum(self.error**2)/len(self.testTrue))
    max_index = list(self.error).index(max(self.error))
    min_index = list(self.error).index(min(self.error))
    max_name = self.test_set[max_index].getName()
    min_name = self.test_set[min_index].getName()
    qtk.report("predicted  MAE", self.MAE)
    qtk.report("predicted RMSE", self.RMSE)
    qtk.report("Maximum  error", max(self.error), max_name)
    qtk.report("Minimum  error", min(self.error), min_name)
    
    def error_estimate(ker, vec):
      tmp = np.vstack([ker, vec])
      K = np.vstack([tmp.T, np.append(vec,1)]).T
      old = np.linalg.det(ker)
      new = np.linalg.det(K)
      
      #kvec = np.dot(ker, vec.T)
      #nvec = np.linalg.norm(vec)
      #nkvec = np.linalg.norm(kvec)
      # angle
      #return np.arccos(np.dot(kvec,vec)/nvec/nkvec)
      # length
      return new/old

    #ee = error_estimate(self.kernelMatrix, self.kernelVectors[0])
    def error_i(i):
      return error_estimate(self.kernelMatrix,
                            self.kernelVectors[i])
    verror = np.vectorize(error_i)
    self.errorEstimate = verror(range(trows))
  ###### END OF PREDICT ######

  def crossValidate():
    pass    

  def singlePrediction(self, i):
    pass
    
