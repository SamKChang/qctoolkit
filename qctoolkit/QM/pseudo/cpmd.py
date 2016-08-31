import qctoolkit as qtk
import re, os
import urllib2
import numpy as np
from qctoolkit.QM.general_io import InpContent

xc_dict = {
            'pbe': 1134,
            'lda': 900,
            'blyp': 1312,
            'factor': 0.6666666667,
            'bp': 1111,
          }

def read(self, path):
  if self.setting['type'] == 'goedecker':
    rev_dict = {v: k for k, v in xc_dict.items()}
    atom = re.compile('^ *&ATOM.*$')
    pot = re.compile('^ *&POTENTIAL.*$')
    end = re.compile('^ *&END.*$')
    try:
      pp_file = urllib2.urlopen(path).readlines()
      pattern = re.compile(r'^.*</*pre>.*$')
      pp_se = filter(pattern.match, pp_file)
      pp_start = pp_file.index(pp_se[0])
      pp_end = pp_file.index(pp_se[1])
      pp_file = pp_file[pp_start:pp_end]
      pp_file[0] = pp_file[0].split('>')[-1]
      for ppStr in pp_file:
        ppStr.replace('&amp;', '&')
    except:
      pp_file = open(path)
    pp = pp_file.readlines()
    i = 0
    while i < len(pp):
      # read ATOM section
      if atom.match(pp[i]):
        while not end.match(pp[i]):
          line = pp[i]
          if re.match('^ *Z *=.*$', line):
            line = re.sub('.*=', '', line)
            self.param['Z'] = float(line.split()[-1])
          elif re.match('^ *ZV *=.*$', line):
            line = re.sub('.*=', '', line)
            self.param['ZV'] = float(line.split()[-1])
          elif re.match('^ *XC *=.*$', line):
            line = re.sub('.*=', '', line)
            key = int(line.split()[0])
            self.param['xc'] = rev_dict[key]
          i += 1
      # read POTENTIAL section
      elif pot.match(pp[i]):
        self.param['l_max'] = int(pp[i+2].split()[0])
        self.param['r_loc'] = float(pp[i+3].split()[0])
        self.param['Cn'] = int(pp[i+4].split()[0])
        Ci = pp[i+4].split()[1:self.param['Cn']+1]
        self.param['Ci'] = [float(c) for c in Ci]
        i += 5
        dim = 0
        self.param['r_nl'] = []
        self.param['h_ij'] = []
        while not end.match(pp[i]):
          line = pp[i]
          data = line.split()
          self.param['r_nl'].append(float(data[0]))
          dim = int(data[1])
          h_ij = np.zeros((dim, dim))
          itr = 2
          for t in range(dim):
            for u in range(t, dim):
              a = t*dim + u
              h_ij[t, u] = float(data[itr])
              h_ij[u, t] = float(data[itr])
              itr += 1
          self.param['h_ij'].append(h_ij)
          i += 1
      i += 1

def write(self, name=None):
  out = InpContent(name)
  out.write('&ATOM\n')
  out.write(' Z  = %4.2f\n' % self.param['Z'])
  out.write(' ZV = %4.2f\n' % self.param['ZV'])
  out.write(' XC = %04d    %12.10f\n'\
             % (xc_dict[self.param['xc']], xc_dict['factor']))
  out.write(' TYPE = NORMCONSERVING GOEDECKER\n')
  out.write('&END\n')
  out.write('&INFO\n')
  out.write(' %s\n' % self.info)
  out.write('&END\n')
  out.write('&POTENTIAL\n')
  out.write('    GOEDECKER\n')
  out.write('  %-33d LMAX\n' % self.param['l_max'])
  out.write('   %12.9f%25s\n' % (self.param['r_loc'], 'RC'))
  out.write('  %d ' % self.param['Cn'])
  for c in self.param['Ci']:
    out.write(' %12.9f ' % c)
  for i in range(len(self.param['h_ij'])):
    h_ij = self.param['h_ij'][i]
    r_nl = self.param['r_nl'][i]
    elem = []
    for j in range(len(h_ij)):
      for k in range(j, len(h_ij)):
        elem.append(h_ij[k, j])
    #upper = [list(v) for v in np.triu(h_ij)]
    #elem = [e for v in upper for e in v if e]
    out.write('\n% 15.9f %2d' % (r_nl, len(h_ij)))
    for e in elem:
      out.write(' %12.9f' % e)
  out.write('\n&END')

  out.close()
