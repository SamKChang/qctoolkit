import qctoolkit as qtk
import numpy as np
import copy
import os
import urllib2
import bigdft
import cpmd
import espresso

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
          'l_max': 0,
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
      if len(path) <= 2:
        self.get(path)
      else:
        self.read(path)
    # sizes of projector matrix
    self.dim = self.getDim()

  def getDim(self):
    n_size = np.array([h.shape[0] for h in self['h_ij']])
    n_dim = 1 + self['Cn'] + len(self['h_ij']) # r_loc Cn and r_nl 
    for i in n_size:
      n = i * (i + 1) / 2
      n_dim += n
    return self['Cn'], len(self['h_ij']), n_dim

  def __repr__(self):
    d = self.dim
    return 'Cn=%d, l_max=%d, parameters=%d\n' % (d[0], d[1], d[2])

  def __getitem__(self, key):
    return self.param[key]

  def __setitem__(self, key, value):
    self.param[key] = value

  def keys(self):
    return self.param.keys()

  def get(self, element, **kwargs):

    def PPName(**kwargs):
      element = kwargs['element'].title()
      if 'vdw' in kwargs:
        PPvdw = '_dcacp_'
      else:
        PPvdw = ''
      if 'd_shell' in kwargs:
        nve = qtk.n2ve(element) + 10
      else:
        nve = qtk.n2ve(element)
      PPStr = element + PPvdw + '_q%d_' % nve +\
        self.setting['xc'].lower() + '.psp'
      return PPStr, element

    def PPCheck(xc, element, pp_file_str, **kwargs):
      if xc == 'lda':
        xc = 'pade'
      ne = qtk.n2ve(element)
      try:
        if 'dcacp' in kwargs and kwargs['dcacp']:
          pp_path = os.path.join(xc.upper(), "%s_DCACP_%s" %\
                    (element, xc.upper()))
          if element in qtk.setting.dcacp_dict:
            pp_path = pp_path + "_%s" % qtk.setting.dcacp_dict[element]
          pp_file = os.path.join(qtk.setting.cpmd_dcacp_url, pp_path)
        else:
          pp_path = os.path.join(xc,
            element + '-q' + str(qtk.n2ve(element)))
          pp_file = os.path.join(qtk.setting.cpmd_pp_url, pp_path)
        saved_pp_path = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
        if not os.path.exists(saved_pp_path) \
        and qtk.setting.download_pp:
          if pp_file:
            new_pp = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
            pp_content = urllib2.urlopen(pp_file).read()
            qtk.report('PPCheck', 'pp file %s not found in %s, ' \
                       % (pp_file_str, qtk.setting.cpmd_pp) + \
                       'but found on internet, download now...')
            new_pp_file = open(new_pp, 'w')
            new_pp_file.write(pp_content)
            new_pp_file.close()
            pp_file = new_pp
        return saved_pp_path
      except:
        qtk.warning('something wrong with pseudopotential')

    kwargs['element'] = element
    file_name, element = PPName(**kwargs)
    if 'vdw' in kwargs and kwargs['vdw'].lower() == 'dcacp':
      dcacp_flag = True
    else:
      dcacp_flag = False
    pp_file = PPCheck(self.setting['xc'], 
              element, 
              file_name, 
              dcacp=dcacp_flag)
    return self.read(pp_file)

  def read(self, path):
    fullpath = os.path.abspath(path)
    self.path, self.name = os.path.split(fullpath)
    if self.setting['program'] == 'cpmd':
      cpmd.read(self, path)
    elif self.setting['program'] == 'bigdft':
      bigdft.read(self, path)

    else:
      qtk.exit('program %s is not implemented for PP'\
        % self.setting['program'])
    return self

  def write(self, name=None, **kwargs):
    if 'inplace' not in kwargs:
      inplace = True
    else:
      inplace = kwargs['inplace']
    if self.setting['program'] == 'cpmd':
      if name:
        stem, ext = os.path.splitext(name)
        if ext != '.psp':
          name = name + '.psp'
      if name and not inplace:
        name = os.path.join(qtk.setting.cpmd_pp, name)
      cpmd.write(self, name)
    elif self.setting['program'] == 'bigdft':
      if name and not inplace:
        name = os.path.join(qtk.setting.bigdft_pp, name)
      bigdft.write(self, name)
    elif self.setting['program'] == 'espresso':
      if name:
        stem, ext = os.path.splitext(name)
        if ext == '.psp' or ext == '.UPF':
          name = stem
      if name and not inplace:
        espresso_name = os.path.join(qtk.setting.espresso_pp, 
          name + '.UPF')
        cpmd_name = os.path.join(qtk.setting.cpmd_pp, 
          name + '.psp')
      else:
        espresso_name = name + '.UPF'
        cpmd_name = name + '.psp'
      espresso.write(self, cpmd_name, espresso_name)

    else:
      qtk.exit('program %s is not implemented for PP'\
        % self.setting['program'])

  def getSize(self):
    return self.param['Cn'], [len(h) for h in self.param['h_ij']]

  def resize(self, cn, hn):
    too_small = False
    if cn < self.param['Cn']:
      too_small = True
    else:
      self.param['Cn'] = cn
      tmp = [0 for i in range(cn)]
      tmp[:len(self.param['Ci'])] = self.param['Ci']
      self.param['Ci'] = copy.deepcopy(tmp)
    if len(hn) < len(self.param['h_ij']):
      too_small = True
    else:
      n = len(hn)
      tmp = [0 for i in range(n)]
      tmp[:len(self.param['r_nl'])] = self.param['r_nl']
      self.param['r_nl'] = copy.deepcopy(tmp)
      tmp = [np.zeros((0,0)) for i in range(len(hn))]
      tmp[:len(self.param['h_ij'])] = self.param['h_ij']
      self.param['h_ij'] = copy.deepcopy(tmp)
    if not too_small:
      for i in range(len(hn)):
        n = hn[i]
        h = self.param['h_ij'][i]
        if n > len(h):
          tmp = np.zeros((n,n))
          tmp[:len(h), :len(h)] = h
          self.param['h_ij'][i] = tmp
    else:
      qtk.warning('PP dimension too small')
    self.dim = self.getDim()
    return self

  def vectorize(self):
    out = [self['ZV']]
    out.append(self['r_loc'])
    for c in self['Ci']:
      out.append(c)
    for i in range(len(self['r_nl'])):
      out.append(self['r_nl'][i])
      for row in range(len(self['h_ij'][i])):
        for col in range(row, len(self['h_ij'][i])):
          out.append(self['h_ij'][i][row, col])
    return np.array(out), self['Cn'], [len(h) for h in self['h_ij']]
       
  def unvectorize(self, param, Cn, h_ij_length, **kwargs):

    def nonZero(value):
      if abs(value) > 1E-8:
        return value
      else:
        return 0.

    self['ZV'] = nonZero(param[0])
    self['r_loc'] = nonZero(param[1])
    self['Ci'] = []
    for i in range(Cn):
      self['Ci'].append(nonZero(param[2+i]))
    self['h_ij'] = []
    itr = 2 + Cn
    self['r_nl'] = []
    for size in h_ij_length:
      self['r_nl'].append(nonZero(param[itr]))
      itr += 1
      h_ij = []
      for i in range(size):
        h_i = []
        for j in range(i):
          h_i.append(0)
        for j in range(i, size):
          h_i.append(nonZero(param[itr]))
          itr += 1
        h_ij.append(h_i)
      h_ij = np.array(h_ij)
      diag = np.diag(np.diag(h_ij))
      h_ij = h_ij + h_ij.T - diag
      self['h_ij'].append(h_ij)
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
