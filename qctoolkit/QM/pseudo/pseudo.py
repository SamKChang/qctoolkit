import qctoolkit as qtk
import cpmd
import numpy as np
import copy
import bigdft

pp_setting = {
               'program': 'cpmd',
               'xc': 'pbe',
               'type': 'goedecker',
             }
param = {
          'Z': 0,
          'ZV': 0,
          'xc': 'pbe',
          'Cn': 0,
          'Ci': [],
          'r_loc': 0,
          'r_nl': [],
          'h_ij': [],
        }

class PP(object):
  def __init__(self, path=None, **kwargs):
    self.param = copy.deepcopy(param)
    self.setting = kwargs
    self.info = ''

    for key, value in pp_setting.iteritems():
      if key in kwargs:
        self.setting[key] = kwargs[key]
      else:
        self.setting[key] = value
    if path:
      self.read(path)

  def __repr__(self):
    return str(self.param)

  def read(self, path):
    if self.setting['program'] == 'cpmd':
      cpmd.read(self, path)
    elif self.setting['program'] == 'bigdft':
      bigdft.read(self, path)

    else:
      qtk.exit('program %s is not implemented for PP'\
        % self.setting['program'])

  def write(self, name=None):
    if self.setting['program'] == 'cpmd':
      cpmd.write(self, name)
    elif self.setting['program'] == 'bigdft':
      bigdft.write(self, name)

    else:
      qtk.exit('program %s is not implemented for PP'\
        % self.setting['program'])

  def resize(self, cn, hn):
    if cn > self.param['Cn']:
      self.param['Cn'] = cn
      tmp = [0 for i in range(cn)]
      tmp[:len(self.param['Ci'])] = self.param['Ci']
      self.param['Ci'] = copy.deepcopy(tmp)
    if len(hn) > len(self.param['h_ij']):
      n = len(hn)
      tmp = [0 for i in range(n)]
      tmp[:len(self.param['r_nl'])] = self.param['r_nl']
      self.param['r_nl'] = copy.deepcopy(tmp)
      tmp = [np.zeros((0,0)) for i in range(len(hn))]
      tmp[:len(self.param['h_ij'])] = self.param['h_ij']
      self.param['h_ij'] = copy.deepcopy(tmp)
    for i in range(len(hn)):
      n = hn[i]
      h = self.param['h_ij'][i]
      if n > len(h):
        tmp = np.zeros((n,n))
        tmp[:len(h), :len(h)] = h
        self.param['h_ij'][i] = tmp
    return self

  def __add__(self, other):
    assert type(other) is type(self)
    assert self.param['xc'] == other.param['xc']

    def vecResize(vec, n):
      tmp = [0 for i in range(n)]
      tmp[:len(vec)] = vec
      return copy.deepcopy(tmp)

    cn = max(self.param['Cn'], other.param['Cn'])
    h1 = [len(h) for h in self.param['h_ij']]
    h2 = [len(h) for h in other.param['h_ij']]
    if len(h1) < len(h2):
      h1 = vecResize(h1, len(h2))
    elif len(h2) < len(h1):
      h2 = vecResize(h2, len(h1))
    hn = [max(i) for i in zip(h1, h2)]
    self.resize(cn, hn)
    other.resize(cn, hn)
    out = copy.deepcopy(self)
    out.param['r_loc'] = self.param['r_loc'] + other.param['r_loc']
    for i in range(len(self.param['Ci'])):
      out.param['Ci'][i] = self.param['Ci'][i] + other.param['Ci'][i]
    for i in range(len(self.param['r_nl'])):
      out.param['r_nl'][i] = self.param['r_nl'][i] +\
                             other.param['r_nl'][i]
    for i in range(len(self.param['h_ij'])):
      h1 = self.param['h_ij'][i]
      h2 = other.param['h_ij'][i]
      out.param['h_ij'][i] = h1 + h2

    out.param['ZV'] = self.param['ZV'] + other.param['ZV']

    return out

  def __mul__(self, other):
    assert type(other) is float or type(other) is int
    out = copy.deepcopy(self)
    out.param['r_loc'] = self.param['r_loc'] * other
    for i in range(len(self.param['Ci'])):
      out.param['Ci'][i] = self.param['Ci'][i] * other
    for i in range(len(self.param['r_nl'])):
      out.param['r_nl'][i] = self.param['r_nl'][i] * other
    for i in range(len(self.param['h_ij'])):
      h1 = self.param['h_ij'][i]
      out.param['h_ij'][i] = h1 * other

    out.param['ZV'] = self.param['ZV'] * other

    return out
