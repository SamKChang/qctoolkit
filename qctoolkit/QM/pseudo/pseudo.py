import qctoolkit as qtk
import cpmd
import numpy as np
import copy
import bigdft
import os
import urllib2

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
    n_size = np.array([h.shape[0] for h in self['h_ij']])
    if 1 in n_size:
      n_dim = sum(n_size * (n_size - 1)) + 2 + len(n_size) + self['Cn']
    else:
      n_dim = sum(n_size * (n_size - 1)) + 1 + len(n_size) + self['Cn']
    self.dim = (self['Cn'], self['l_max'], n_dim)

  def __repr__(self):
    d = self.dim
    return 'Cn=%d, l_max=%d, parameters=%d\n' % (d[0], d[1], d[2])

  def __getitem__(self, key):
    return self.param[key]

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

    # not used by PP object but by QMInp cpmd parts
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
            print new_pp
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
