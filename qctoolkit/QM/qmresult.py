import qctoolkit as qtk
import numpy as np
import pandas as pd
import glob, os, re, copy

class QMResult(object):
  def __init__(self, pattern, qm_property='Et', **kwargs):

    if 'program' not in kwargs:
      kwargs['program'] = qtk.setting.qmcode
    if 'unit' not in kwargs:
      kwargs['unit'] = 'Eh'

    out_files = sorted(glob.glob(pattern))
    if not out_files:
      qtk.warning('no output files matched %s' % pattern)

    files = []
    values = []
    self.qmout = []
    self.names = []
    for out in out_files:
      out_file = os.path.split(out)[1]
      out_file = re.sub('\.out', '', out_file)
      qmout = qtk.QMOut(out, program=kwargs['program'])
      qmout.inUnit(kwargs['unit'])
      try:
        value = getattr(qmout, qm_property)
      except:
        qtk.warning('%s not found for file %s' %\
          (qm_property, out_file))
        value = np.nan
      files.append(out_file)
      values.append(value)
      self.names.append(out_file)
      self.qmout.append(qmout)
    self.data = pd.Series(values, index=files)

    method_strlist = [m for m in dir(self.data)]
    method_list = []
    for m in method_strlist:
      try:
        if callable(getattr(self.data, m)):
          method_list.append(m)
      except:
        pass

    p = re.compile('^[a-z].*$')
    method_list = filter(p.match, method_list)
    for m in method_list:
      if m not in dir(self):
        setattr(self, m, getattr(self.data, m))

  def __repr__(self):
    return str(self.data)

  def __add__(self, other):
    out = copy.deepcopy(self)
    if type(other) is type(self):
      out.data = self.data + other.data
    elif type(other) is int or type(other) is float:
      out.data = self.data + other
    return out

  def __sub__(self, other):
    out = copy.deepcopy(self)
    if type(other) is type(self):
      out.data = self.data - other.data
    elif type(other) is int or type(other) is float:
      out.data = self.data - other
    return out

  def __mul__(self, other):
    assert type(other) is int or type(other) is float
    out = copy.deepcopy(self)
    if type(other) is int or type(other) is float:
      out.data = self.data * other
    return out

  def __div__(self, other):
    assert type(other) is int or type(other) is float
    out = copy.deepcopy(self)
    if type(other) is int or type(other) is float:
      out.data = self.data / other
    return out

  def indexAdd(self, factor):
    index_list = self.data.index.tolist()
    index = [float(i)+factor for i in index_list]
    new_index = pd.Index(index)
    self.data.index = new_index

  def indexSub(self, factor):
    index_list = self.data.index.tolist()
    index = [float(i)-factor for i in index_list]
    new_index = pd.Index(index)
    self.data.index = new_index

  def indexMul(self, factor):
    index_list = self.data.index.tolist()
    index = [float(i)*factor for i in index_list]
    new_index = pd.Index(index)
    self.data.index = new_index

  def indexDiv(self, factor):
    index_list = self.data.index.tolist()
    index = [float(i)/factor for i in index_list]
    new_index = pd.Index(index)
    self.data.index = new_index

  def indexExtract(self, pattern):
    index_list = self.data.index.tolist()
    p = re.compile(pattern)
    index = [p.match(i).group(1) for i in index_list]
    new_index = pd.Index(index)
    self.data.index = new_index

  def indexRemove(self, pattern, **kwargs):
    index_list = self.data.index.tolist()
    index = [re.sub(pattern, '', i) for i in index_list]
    new_index = pd.Index(index)
    self.data.index = new_index
