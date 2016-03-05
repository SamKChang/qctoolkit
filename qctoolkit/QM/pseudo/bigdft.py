import qctoolkit as qtk
import re, os
import urllib2
import datetime
import numpy as np
from qctoolkit.QM.general_io import InpContent

type_dict = {
              'goedecker': 10,
              'nlcc': 12,
            }
xc_dict = {
           'ldaxc': '-020',
           'pbexc': '-101130',
           'blypxc': '-106131',
           'lda': '1',
           'pbe': '11',
           'blyp': '18',
           'bp': '19',
         }

def read(self, path):
  try:
    pp_file = urllib2.urlopen(path).readlines()
  except:
    pp_file = open(path)
  pp = pp_file.readlines()
  rev_xc = {v: k for k, v in xc_dict.items()}
  rev_type = {v: k for k, v in type_dict.items()}
  self.param['Z'] = float(pp[1].split()[0])
  self.param['ZV'] = float(pp[1].split()[1])
  self.setting['type'] = rev_type[int(pp[2].split()[0])]
  self.param['xc'] = rev_xc[pp[2].split()[1]]
  rc_line = pp[3].split()
  self.param['r_loc'] = float(rc_line[0])
  self.param['Cn'] = int(rc_line[1])
  for i in range(self.param['Cn']):
    self.param['Ci'].append(float(rc_line[i+2]))
  self.param['l_max'] = int(pp[4].split()[0])

  shift = 5
  i = 0
  while i < self.param['l_max']:
  #for i in range(self.param['l_max']):
    line = pp[shift].split()
    self.param['r_nl'].append(float(line[0]))
    size = int(float(line[1]))
    h_data = line[2:]
    h_jk = np.zeros((size, size))
    for j in range(size):
      for k in range(j, size):
        h_jk[j,k] = h_data[k-j]
      shift = shift + 1
      h_data = pp[i+shift].split()
    self.param['h_ij'].append(h_jk)
    i = i + 1
  if self.setting['type'] == 'nlcc':
    shift = shift  + len(h_jk)
    self.param['rcore'] = float(pp[shift].split()[0])
    self.param['qcore'] = float(pp[shift].split()[1])

def write(self, name=None):
  today = datetime.date.today()
  out = InpContent(name)
  out.write('Goedecker pseudopotential for %s\n' % \
    qtk.Z2n(self.param['Z']))
  out.write('   %2d  %2d  %s ! Z, ZV, date(ddmmyy)\n' %\
    (self.param['Z'], self.param['ZV'], today.strftime('%d%m%y')))
  out.write('   %d  %s  2 0 2002 0  ! xc:%s\n' %\
    (type_dict[self.setting['type']], 
     xc_dict[self.param['xc']], 
     self.param['xc']))
  out.write('   %12.8f %3d' % (self.param['r_loc'], self.param['Cn']))
  for i in range(self.param['Cn']):
    out.write(' % 12.8f' % self.param['Ci'][i])
  out.write(' ! rloc, #C, C[i]')
  out.write('\n  %d\n' % self.param['l_max'])
  for i in range(len(self.param['r_nl'])):
    r_nl = self.param['r_nl'][i]
    h_ij = self.param['h_ij'][i]
    out.write('   %12.8f %3d' % (r_nl, len(h_ij)))
    if len(h_ij) > 0:
      for j in range(len(h_ij)):
        if j > 0:
          out.write('%19s' % '')
        for k in range(j, len(h_ij)):
          out.write(' % 12.8f' % h_ij[j][k])
        if j == 0:
          out.write(' ! h_ij\n')
        else:
          out.write('\n')
    else:
      out.write(' % 12.8f\n' % 0)
  if self.param['r_nl'] > 0:
    for j in range(len(h_ij)):
      out.write('%19s' % '')
      for k in range(j, len(h_ij)):
        out.write(' % 12.8f' % 0.0)
      if j == 0:
        out.write(' ! abinit k_ij\n')
      else:
        out.write('\n')
  if self.setting['type'] == 'nlcc':
    out.write('   %12.8f     %12.8f' %\
      (self.param['rcore'], self.param['qcore']))
  out.close()
