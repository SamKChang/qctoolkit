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
        'IOp(2/12=3)', # allow atoms to be too near
        'force',
      ]
      self.setting['gaussian_setting'] = gaussian_setting

  def run(self, name=None, **kwargs):
    self.setting.update(kwargs)
    return univ.runCode(self, GaussianBasisInput, name, **self.setting)
    
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
    if 'cc' in theory:
      try:
        ind = self.setting['gaussian_setting'].index('force')
        del self.setting['gaussian_setting'][ind]
      except:
        pass
      
    if 'openshell' in self.setting :
      if self.setting['openshell'] == 'restricted':
        theory = 'ro' + theory
      elif self.setting['openshell'] == 'unrestricted':
        theory = 'u' + theory
    if 'geopt' in self.setting and self.setting['geopt']:
      self.setting['gaussian_setting'].append('opt')
      if 'force' in self.setting['gaussian_setting']:
        ind = self.setting['gaussian_setting'].index('force')
        del self.setting['gaussian_setting'][ind]
    if 'print_polarizability' in self.setting\
    and self.setting['print_polarizability']:
      self.setting['gaussian_setting'].append('polar')
      if 'force' in self.setting['gaussian_setting']:
        ind = self.setting['gaussian_setting'].index('force')
        del self.setting['gaussian_setting'][ind]
    basis = self.setting['basis_set']
    if type(basis) is str:
      if 'def2' in basis.lower():
        basis = basis.replace('-', '')
    if 'charge_multiplicity' not in self.setting:
      charge, multiplicity = \
        self.molecule.charge, self.molecule.multiplicity
    else:
      charge, multiplicity = \
        self.setting['charge_multiplicity']
      if (charge + sum(self.molecule.Z) + 1) % 2 != multiplicity % 2:
        qtk.exit("charge: %3.1f and multiplicity %d is " \
                 % (float(charge), multiplicity) + 
                 "not competible with number of electron %d\n"
                 % sum(self.molecule.Z))

    gaussian_setting = self.setting['gaussian_setting']

    chk_flag = False
    save_list = [
      'save_density',
      'ks_states',
      'save_wf',
      'save_restart',
    ]
    for s in save_list:
      if s in self.setting and self.setting[s]:
        chk_flag = True

    if 'threads' in self.setting:
      inp.write('%%nproc=%d\n' % self.setting['threads'])
    else:
      inp.write('%nproc=\n')
    if chk_flag:
      if name:
        inp.write('%%chk=%s.chk\n' % name)
      else:
        inp.write('%chk=tmp.chk\n')

      density_dict = {'ccsd', 'mp2', 'mp3', 'ccd', 'cid', 'cisd'}

      if 'save_density' in self.setting\
      and self.setting['save_density']:
        if theory.lower() in density_dict\
        and 'Density=Current' not in self.setting['gaussian_setting']:
          self.setting['gaussian_setting'].append('Density=Current')
       
    if 'nuclear_charges' in self.setting:
      gaussian_setting.append('Charge')
    if 'vdw' in self.setting:
      if self.setting['vdw'] == 'd3':
        gaussian_setting.append('EmpiricalDispersion=GD3')
      elif self.setting['vdw'] == 'd2':
        gaussian_setting.append('EmpiricalDispersion=GD2')
    if 'print_energy' in self.setting and self.setting['print_energy']:
      gaussian_setting.append('ExtraLinks=L608')
    if 'mode' in self.setting and self.setting['mode'] == 'geopt':
      gaussian_setting.append('Opt')
    if type(basis) is str:
      inp.write("# %s/%s" % (theory, basis))
    else:
      inp.write("# %s/gen" % theory)
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

    if type(basis) is dict:
      inp.write('\n')
      l_dict = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G'}
      for atom in set(molecule.type_list):
        if atom not in basis:
          qtk.exit("basis function for atom %s is missing" % atom)
        inp.write("%-2s   0\n" % atom)
        for b in basis[atom]:
          inp.write("%s  %2d   1.00\n" % (l_dict[b[0]], len(b)-1))
          for ec in b[1:]:
            inp.write("  %16.8f  %16.8f\n" % (ec[0], ec[1]))
        inp.write('****\n')
    inp.write('\n\n')
    inp.close()

    return inp

class out(GaussianBasisOutput):
  def __init__(self, qmout=None, **kwargs):
    GaussianBasisOutput.__init__(self, qmout, **kwargs)
    if qmout:
      outfile = open(qmout)
      data = outfile.readlines()

      pattern = re.compile(" *R *M *S *D *=")
      try:
        report = filter(lambda x: "\\" in x, data)
      except Exception as e:
        qtk.exit("error when access final report with message %s" % str(e))
      final_str = ''.join(report)
      final_str = final_str.replace('\n', '')
      final_list = final_str.split('\\')
      try:
        rmsd = filter(pattern.match, final_list)[0]
      except Exception as e:
        qtk.exit("something wrong when accessing final energy" +\
          " with error message: %s" % str(e))
        
      ind = final_list.index(rmsd) - 1
      Et_str = final_list[ind]
      term = Et_str.split('=')[0].replace(' ', '')
      E_list = ['HF', 'CCSD', 'CCSD(T)', 'MP4SDQ', 'MP4DQ'
                'MP4D', 'MP3', 'MP2', 'HF']
      while term not in E_list:
        ind = ind - 1
        Et_str = final_list[ind]
        term = Et_str.split('=')[0].replace(' ', '')
        
      self.Et = float(Et_str.split('=')[1].replace(' ',''))
      self.detail = final_list

      EJStr = filter(lambda x: 'EJ=' in x, data)
      if len(EJStr) > 0:
        EJStr = EJStr[-1]
        EJind = data.index(EJStr)
        EJStr = data[EJind].replace(' ', '')
        EComponents = filter(None, EJStr.split('E'))
        tags_dict = {
          'T':'Ekin',
          'V':'Eext',
          'J':'Eee',
          'K':'Ex',
          'Nuc':'Enn',
        }
        for EStr in EComponents:
          tag, E = EStr.split('=')
          try:
	          self.energies[tags_dict[tag]] = float(E)
          except Exception as e:
            qtk.warning("FORTRAN float overflow " + \
                        "when extracting energy components for " +\
                        qmout + " with error message %s" % str(e))
            self.energies[tags_dict[tag]] = np.nan

      crdStr = filter(lambda x: 'Angstroms' in x, data)[-1]
      ind = len(data) - data[::-1].index(crdStr) + 2
      ZR = []
      while True:
        if not data[ind].startswith(' ---'):
          crdData = [float(c) 
                     for c in filter(None, data[ind].split(' '))]
          crd = [crdData[1]]
          crd.extend(crdData[3:])
          ZR.append(crd)
          ind = ind + 1
        else:
          break
      self.molecule = qtk.Molecule()
      self.molecule.build(ZR)
      self.nuclear_repulsion = self.molecule.nuclear_repulsion()

      force = []
      fStr_list = filter(lambda x: 'Forces (Hartrees' in x, data)
      if len(fStr_list) > 0:
        fStr = fStr_list[-1]
        ind = len(data) - data[::-1].index(fStr) + 2
        for i in range(self.molecule.N):
          fLst = filter(None, data[ind+i].split(' '))[2:]
          force.append([float(s) for s in fLst])
        self.force = np.array(force)
      else:
        self.force = np.nan

      dipole = []
      uStr_list = filter(lambda x: 'Debye)' in x, data)
      if len(uStr_list) > 0:
        uStr = uStr_list[-1]
        ind = len(data) - data[::-1].index(uStr)
        dipoleStr = filter(None, data[ind].split(' '))
        for i in [1,3,5]:
          dipole.append(float(dipoleStr[i]))
        self.dipole = np.array(dipole)
      else:
        self.dipole = np.nan

      qp = []
      qStr_list = filter(lambda x: ' Quadrupole moment' in x, data)
      if len(qStr_list) > 0:
        qStr = qStr_list[-1]
        ind = len(data) - data[::-1].index(qStr)
        for i in range(2):
          tmp = dipoleStr = filter(None, data[ind+i].split(' '))
          for j in [1,3,5]:
            qp.append(float(tmp[j]))
        xx, yy, zz, xy, xz, yz = qp
        self.quadrupole = np.array([
          [xx, xy, xz],
          [xy, yy, yz],
          [xz, yz, zz],
        ])
      else:
        self.quadrupole = np.nan

      read_fchk = True
      if 'read_fchk' in kwargs:
        read_fchk = kwargs['read_fchk']

      if read_fchk:
        fchk = os.path.join(self.path, self.stem) + ".fchk"
        if os.path.exists(fchk):
          if 'debug' in kwargs and kwargs['debug']: 
            self.getMO(fchk)
          else:
            try:
              self.getMO(fchk)
            except Exception as e:
              qtk.warning("something wrong while loading fchk file"+\
                " with error message: %s" % str(e))

  def getMO(self, fchk):
    self.program = 'gaussian'
    fchkfile = open(fchk)
    fchk = fchkfile.readlines()
    self.basis_name = filter(None, fchk[1].split(' '))[2]

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
    try:
      self.mo_eigenvalues_beta, _ = \
        readFchk('Beta Orbital Energies')
      _mob, _dimb = readFchk('Beta MO coefficients')
    except:
      _mob, _dimb = None, None
      self.mo_eigenvalues_beta = None
    self.n_ao = _dim/self.n_mo
    self.n_basis = self.n_ao
    self.mo_vectors = np.array(_mo).reshape([self.n_mo, self.n_ao])
    if self.mo_eigenvalues_beta is not None:
      self.mo_vectors_beta = np.array(_mob).reshape([self.n_mo, self.n_ao])

    _map = readFchk('Shell to atom map', int)[0]
    _map_coord = list(np.diff(_map))
    _map_coord.insert(0,1)
    _R_ind = [i for i in range(_n_shell) if _map_coord[i]>0]
    self.R_bohr = np.array([_coord[i] for i in _R_ind])
    self.R = self.R_bohr / 1.88972613
    _neStr = filter(lambda x: 'Number of electrons' in x, fchk)[0]
    _ne = float(filter(None, _neStr.split(' '))[-1])
    self.occupation = []
    warned = False
    for i in range(self.n_mo):
      if i < _ne/2 - 1:
        self.occupation.append(2.0)
      elif i >= _ne/2 -1 and i < _ne/2:
        if abs(_ne % 2) < 1E-5:
          self.occupation.append(2.0)
        else:
          self.occupation.append(_ne % 2)
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

      if _types[i] == 0:
        l_max = 1
        blist = ['s']
      else:
        if _types[i] > 0:
          l_max = sum(range(_types[i]+2))
          blist = basis_list[_types[i]]
        else:
          l_max = 2 * abs(_types[i]) + 1
          blist = ['D0' for _ in range(l_max)]
          if not warned:
            qtk.warning("sp shell or spherical basis are used, " +\
              "some functions will not work")
            warned = True
 
      for l in range(l_max):
        bfn = copy.deepcopy(bfn_base)
        bfn['exponents'] = copy.deepcopy(exp)
        bfn['coefficients'] = copy.deepcopy(cef)
        bfn['type'] = blist[l]
        tmp = tmp + 1
        self.basis.append(bfn)
    self.basisFormat()
