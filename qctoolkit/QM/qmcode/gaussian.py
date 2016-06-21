import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import periodictable as pt
import universal as univ

class inp(GaussianBasisInput):
  def __init__(self, molecule, **kwargs):
    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()
    if 'gaussian_setting' not in self.setting:
      gaussian_setting = [
        '6d 10f',
        'nosymm',
        'Scf(maxcycle=1000,verytight)',
        'int(grid=ultrafine)',
      ]
      self.setting['gaussian_setting'] = gaussian_setting

  def run(self, name=None, **kwargs):
    pass
#    # creat_folder routine check default name
#    name = self.create_folder(name)
#    cwd = os.getcwd()
#    os.chdir(name)
#    inp = name + '.inp'
#    self.write(inp)
#    try:
#      out = qmjob.QMRun(inp, 'gaussian', **kwargs)
#    except:
#      qtk.warning("qmjob finished unexpectedly for '" + \
#                  name + "'")
#      out = GaussianBasisOutput(program='gaussian')
#    finally:
#      os.chdir(cwd)
#    return out
    
  def write(self, name=None, **kwargs):
    self.setting.update(kwargs)
    univ.cornerCube(self)
    inp, molecule = \
      super(GaussianBasisInput, self).write(name, **self.setting)

    theory_dict = {
      'pbe': 'pbepbe',
      'pbe0': 'pbe1pbe',
    }

    if self.setting['theory'] in theory_dict:
      theory = theory_dict[self.setting['theory']]
    else:
      theory = self.setting['theory']
    basis = self.setting['basis_set']
    if 'def2' in basis.lower():
      basis = basis.replace('-', '')
    charge, multiplicity = \
      self.molecule.charge, self.molecule.multiplicity

    gaussian_setting = self.setting['gaussian_setting']

    chk_flag = False
    save_list = [
      'save_density',
      'ks_states',
      'save_wf',
    ]
    for s in save_list:
      if s in self.setting and self.setting[s]:
        chk_flag = True

    if 'threads' in self.setting:
      inp.write('%%nproc=%d\n' % self.setting['threads'])
    else:
      inp.write('%nproc=\n')
    if chk_flag:
      inp.write('%chk=\n')

      density_dict = {'ccsd', 'mp2', 'mp3', 'ccd', 'cid', 'cisd'}

      if 'save_density' in self.setting\
      and self.setting['save_density']\
      and theory.lower() in density_dict\
      and 'Density=Current' not in self.setting['gaussian_setting']:
        self.setting['gaussian_setting'].append('Density=Current')
    if 'nuclear_charges' in self.setting:
      gaussian_setting.append('Charge')
    inp.write("# %s/%s" % (theory, basis))
    for s in list(set(gaussian_setting)):
      inp.write(" %s" % s)
    inp.write("\n\n%s\n\n" % self.molecule.name)
    inp.write("%d   %d\n" % (charge, multiplicity))

    for i in range(molecule.N):
      inp.write('%-2s % 8.4f % 8.4f % 8.4f\n' % \
                (molecule.type_list[i],
                 molecule.R[i, 0],
                 molecule.R[i, 1],
                 molecule.R[i, 2]))

    if 'nuclear_charges' in self.setting:
      new_Z = self.setting['nuclear_charges']
    else:
      new_Z = []
    if new_Z:
      inp.write('\n')
      for chg_list in new_Z:
        for i in range(molecule.N):
          if chg_list[0] == i:
            Ri = molecule.R[i]
            charge = chg_list[1]
            inp.write('% 8.4f % 8.4f % 8.4f  % .3f\n' %\
                      (Ri[0], Ri[1], Ri[2], charge))
    inp.write('\n\n')
    inp.close()

    return inp

class out(GaussianBasisOutput):
  def __init__(self, qmout=None, **kwargs):
    GaussianBasisOutput.__init__(self, qmout, **kwargs)
    if qmout:
      outfile = open(qmout)
      data = outfile.readlines()

      report_str = filter(lambda x: 'N-N=' in x, data)[0]
      start = data.index(report_str)
      not_found = True
      i = 1
      while not_found:
        end = start + i
        i += 1
        if data[end] == '\n':
          not_found = False
          break
      report = data[start:end]
      final_str = ''.join(report)
      final_str = final_str.replace('\n', '')
      final_list = final_str.split('\\')
      pattern = re.compile("R *M *S *D *=")
      rmsd = filter(pattern.match, final_list)[0]
      ind = final_list.index(rmsd) - 1
      Et_str = final_list[ind]
      self.Et = float(Et_str.split('=')[1].replace(' ',''))
      self.detail = final_list
      
      fchk = os.path.join(self.path, self.stem) + ".fchk"
      if os.path.exists(fchk):
        if qtk.setting.debug: 
          self.getMO(fchk)
        else:
          try:
            self.getMO(fchk)
          except:
            qtk.warning("something wrong while loading fchk file")

  def getMO(self, fchk):
    fchkfile = open(fchk)
    fchk = fchkfile.readlines()

    def basisList(L):
      N_dict = {4: 'g', 5: 'h', 6: 'k'}
      orbital = ['x', 'y', 'z']
      base = [2 for _ in range(L)]

      out = []
      for n in range(L + 1):
        b = copy.deepcopy(base)
        for i in range(n):
          b[i] = 0
        orb = N_dict[L] + ''.join([orbital[j] for j in b])
        out.append(orb)
        for i in range(n, L):
          b[i] = 1
          orb = N_dict[L] + ''.join([orbital[j] for j in b])
          out.append(orb)

      return out

    basis_list = [
      ['s'],
      ['px', 'py', 'pz'],
      ['dxx', 'dyy', 'dzz', 'dxy', 'dxz', 'dyz'],
      ['fxxx', 'fyyy', 'fzzz', 'fxyy', 'fxxy', 
       'fxxz', 'fxzz', 'fyzz', 'fyyz', 'fxyz'],
      basisList(4),
      basisList(5),
    ]

    def readFchk(flag, type=float):
      flagStr = filter(lambda x: flag in x, fchk)[0]
      n_entry = int(filter(None, flagStr.split(' '))[-1])
      ind = fchk.index(flagStr) + 1
      if type is float:
        factor = 5
      elif type is int:
        factor = 6
      n_lines = n_entry / factor + 1
      if not n_entry % factor: n_lines = n_lines - 1
      data = []
      for i in range(ind, ind + n_lines):
        dList = list(np.array(filter(None, 
                     fchk[i].split(' '))).astype(type))
        data.extend(dList)
      return data, n_entry

    self.Z, self.N = readFchk('Nuclear charges')
    self.type_list = [qtk.Z2n(z) for z in self.Z]
    _types, _n_shell = readFchk('Shell types', int)
    _exp, _nfn = readFchk('Primitive exponents')
    _cef = readFchk('Contraction coefficients')[0]
    _coord, _test= readFchk('Coordinates of each shell')
    _ng = readFchk('Number of primitives per shell', int)[0]
    _coord = np.array(_coord).reshape([_n_shell, 3])
    self.mo_eigenvalues, self.n_mo = \
      readFchk('Alpha Orbital Energies')
    _mo, _dim = readFchk('Alpha MO coefficients')
    _map = readFchk('Shell to atom map', int)[0]
    self.n_ao = _dim/self.n_mo
    self.n_basis = self.n_ao
    self.mo_vectors = np.array(_mo).reshape([self.n_mo, self.n_ao])
    _map_coord = list(np.diff(_map))
    _map_coord.insert(0,1)
    _R_ind = [i for i in range(_n_shell) if _map_coord[i]>0]
    self.R_bohr = np.array([_coord[i] for i in _R_ind])
    self.R = self.R_bohr / 1.88972613
    _neStr = filter(lambda x: 'Number of electrons' in x, fchk)[0]
    _ne = float(filter(None, _neStr.split(' '))[-1])
    self.occupation = []
    for i in range(self.n_mo):
      if i < _ne/2:
        self.occupation.append(2.0)
      else:
        self.occupation.append(0.0)
      

    tmp = 0
    shift = 0
    self.basis = []
    for i in range(_n_shell):
      _ind = _map[i] - 1
      bfn_base = {}
      bfn_base['index'] = i
      bfn_base['center'] = _coord[i]
      bfn_base['atom'] = self.type_list[_ind]

      exp = []
      cef = []
      for g in range(_ng[i]):
        exp.append(_exp[i+shift+g])
        cef.append(_cef[i+shift+g])
      shift = shift + _ng[i] - 1
 
      for l in range(sum(range(_types[i]+2))):
        bfn = copy.deepcopy(bfn_base)
        bfn['exponents'] = copy.deepcopy(exp)
        bfn['coefficients'] = copy.deepcopy(cef)
        bfn['type'] = basis_list[_types[i]][l]
        tmp = tmp + 1
        self.basis.append(bfn)
