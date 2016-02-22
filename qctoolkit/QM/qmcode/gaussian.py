import qctoolkit as qtk
from qctoolkit.QM.gaussianbasis_io import GaussianBasisInput
from qctoolkit.QM.gaussianbasis_io import GaussianBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob
import periodictable as pt

class inp(GaussianBasisInput):
  def __init__(self, molecule, **kwargs):
    GaussianBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)
    self.backup()

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
    
  def write(self, name=None):
    pass
#    self.reset()
#    molecule = copy.deepcopy(self.molecule)
#    self.cm_check(molecule)
#    if name:
#      nameStr = re.sub('\.inp', '', name)
#      out = os.path.splitext(name)
#      if out[1] != '.inp':
#        name = out[0] + '.inp'
#      if os.path.exists(name) and not qtk.setting.no_warning:
#        qtk.prompt(name + ' exists, overwrite?')
#        try:
#          os.remove(name)
#        except:
#          qtk.exit("can not remove file: " + name)
#    else:
#      nameStr = re.sub('\.inp', '', molecule.name)
#    inp = sys.stdout if not name else open(name, 'w')
#
#    # mode setup
#    if self.setting['mode'] == 'single_point':
#      operation = 'energy'
#    elif self.setting['mode'] == 'geopt':
#      operation = 'optimize'
#
#    # Theory setupt
#    dft = {
#      'pbe': 'xpbe96 cpbe96',
#      'pbe0': 'pbe0',
#      'blyp': 'becke88 lyp',
#      'b3lyp': 'b3lyp',
#      'bp91': 'becke88 perdew91',
#      'bp86': 'becke88 perdew86',
#      'pw91': 'xperdew91 perdew91',
#    }
#    scf = {
#      'rhf', 'rohf', 'uhf',
#    }
#    tce = {
#      'mp2', 'ccsd', 'ccsd(t)',
#    }
#    if self.setting['theory'] in dft:
#      module = 'dft'
#      xc = dft[self.setting['theory']]
#    elif self.setting['theory'] in scf:
#      module = 'scf'
#    elif self.setting['theory'] in tce:
#      module = 'tce'
#
#    # print file
#    inp.write('title %s\n' % nameStr)
#    inp.write('start %s\n' % nameStr)
#    inp.write('echo\n')
#    if self.setting['fix_molecule']\
#    and self.setting['fix_molecule']:
#      fix_mol = 'nocenter'
#    else:
#      fix_mol = ''
#    inp.write('\ngeometry units angstrom %s\n' % fix_mol)
#    if fix_mol:
#      inp.write('symmetry group c1\n')
#    if 'nuclear_charges' in self.setting:
#      new_Z = self.setting['nuclear_charges']
#    else:
#      new_Z = []
#    z_itr = 0
#    for i in range(molecule.N):
#      n_charge = ''
#      e_str = molecule.type_list[i]
#      if new_Z and z_itr < len(new_Z):
#        if new_Z[z_itr][0] == i+1:
#          n_charge = 'charge %.3f' % float(new_Z[z_itr][1])
#          e_str = new_Z[z_itr][2]
#          molecule.string[i] = e_str
#          z_itr = z_itr + 1
#      inp.write(' %-8s % 8.4f % 8.4f % 8.4f %s\n' % \
#                (e_str,
#                 molecule.R[i, 0],
#                 molecule.R[i, 1],
#                 molecule.R[i, 2],
#                 n_charge
#               ))
#    inp.write('end\n\n')
#    inp.write('basis\n')
#    if self.setting['basis_set'] != 'gen':
#      eStr = []
#      for e in range(len(molecule.Z)):
#        if molecule.string[e]:
#          eString = molecule.string[e]
#        else:
#          eString = molecule.type_list[e]
#        if eString not in set(eStr):
#          eStr.append(eString)
#          inp.write(' %-8s library %2s %s\n' % (\
#            eString,
#            molecule.type_list[e],
#            self.setting['basis_set'].upper()),
#          )
#    inp.write('end\n\n')
#    if module == 'dft':
#      inp.write('dft\n')
#      inp.write(' xc %s\n' % xc)
#      inp.write('end\n\n')
#    elif module == 'scf':
#      pass
#    if molecule.charge != 0:
#      inp.write('charge % 5.2f\n' % molecule.charge)
#    if molecule.multiplicity != 1:
#      inp.write('multiplicity %d\n' % molecule.multiplicity)
#
#    inp.write('\ntask %s %s\n' % (module, operation))
#
#    if name:
#      inp.close()


class out(GaussianBasisOutput):
  def __init__(self, qmout=None, **kwargs):
    GaussianBasisOutput.__init__(self, qmout, **kwargs)
    if qmout:
      outfile = open(qmout)
      data = outfile.readlines()
      fchk = os.path.join(self.path, self.stem) + ".fchk"
      if os.path.exists(fchk):
        self.getMO(fchk)

  def getMO(self, fchk):
    fchkfile = open(fchk)
    fchk = fchkfile.readlines()

    basis_list =[
      ['s'],
      ['px', 'py', 'pz'],
      ['dxx', 'dyy', 'dzz', 'dxy', 'dxz', 'dyz'],
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

#        chargeStr = filter(lambda x: 'Nuclear charges' in x,
#                           fchk)[0]
#        self.N = int(filter(None, chargeStr.split(' '))[-1])
#        chargeInd = fchk.index(chargeStr) + 1
#        cLines = self.N / 5 + 1
#        if not self.N % 5: cLine = cLine - 1
#        self.Z = []
#        for i in range(chargeInd, chargeInd + cLines):
#          _list = list(np.array(filter(None, 
#                    fchk[i].split(' '))).astype(float))
#          self.Z.extend(_list)
#        print self.Z


#  """
#  directly parse vasp xml output, 'vasprun.xml'
#  converged energy, system info, and scf steps are extracted
#  """
#  def __init__(self, qmout=None, **kwargs):
#    GaussianBasisOutput.__init__(self, qmout=None, **kwargs)
#    if qmout:
#      stem, ext = os.path.splitext(qmout)
#      outfile = open(qmout, 'r')
#      data = outfile.readlines()
#      Et = filter(lambda x: 'Total' in x and 'energy' in x, data)
#      try:
#        self.Et = float(Et[-1].split( )[-1])
#      except:
#        self.Et = np.nan
#      # necessary for opening *.movecs file
#      n_basis = filter(lambda x: ' functions ' in x, data)
#      try:
#        self.n_basis = int(n_basis[-1].split('=')[1])
#      except:
#        self.n_basis = np.nan
#
#      ######################################
#      # extract basis function information #
#      ######################################
#      basis_P = re.compile(r"  [0-9] [A-Z] +")
#      batom_P = re.compile(r"^  [A-Za-z\-_\.]* *\([A-Z][a-z]*\)")
#      bname_P = re.compile(r"\((.*)\)")
#      coord_P = re.compile(r"^ [A-Za-z\.\-_]+ +[- ][0-9\.]{9,}" +\
#                           r" +[- ][0-9\. ]+$")
#      basisStr = filter(basis_P.match, data)
#      batomStr = filter(batom_P.match, data)
#      coordStr = filter(coord_P.match, data)
#
#      _exponents = [float(filter(None, s.split(' '))\
#        [2]) for s in basisStr]
#      _coefficients = [float(filter(None, s.split(' '))\
#        [3]) for s in basisStr]
#      _N = [int(filter(None, s.split(' '))[0])\
#        for s in basisStr]
#      _type = [filter(None, s.split(' '))[1]\
#        for s in basisStr]
#      _bfnInd = [data.index(batom) for batom in batomStr]
#      _bfnEndPtn = re.compile(r" Summary of \"")
#      _bfnEndStr = filter(_bfnEndPtn.match, data)[0]
#      _bfnInd.append(data.index(_bfnEndStr))
#
#      _ao_keys = [0]
#      for ind in range(len(_bfnInd)-1):
#        _key = _ao_keys[-1]
#        for i in range(_bfnInd[ind]+4, _bfnInd[ind+1]):
#          if len(data[i]) > 1:
#            _key = _key + 1
#        _ao_keys.append(_key)
#      _atoms = [getattr(pt, bname_P.match(
#        filter(None, s.split(' '))[1]).group(1).lower()).symbol\
#        for s in batomStr]
#      self.type_list = [re.split(r'[\._]',
#        filter(None, s.split(' '))[0])[0].title()\
#        for s in coordStr]
#      self.R = np.array([filter(None, s.split(' '))[1:4]\
#        for s in coordStr]).astype(float)
#      self.N = len(self.R)
#      self.Z = [qtk.n2Z(e) for e in self.type_list]
#      self.R_bohr = 1.889725989 * self.R
#
#      _N.append(0)
#      self.basis = []
#      for i in range(len(self.type_list)):
#        e = self.type_list[i]
#        center = self.R_bohr[i]
#        ind = _atoms.index(e)
#        bfn_base = {}
#        bfn_base['atom'] = e
#        bfn_base['center'] = center
#        bfn_base['index'] = i
#        exp = []
#        cef = []
#        for g in range(_ao_keys[ind], _ao_keys[ind+1]):
#          exp.append(_exponents[g])
#          cef.append(_coefficients[g])
#          if _N[g] != _N[g+1] or g+1 >= _ao_keys[ind+1]:
#            bfn = copy.deepcopy(bfn_base)
#            bfn['exponents'] = copy.deepcopy(exp)
#            bfn['coefficients'] = copy.deepcopy(cef)
#            if _type[g] == 'S':
#              bfn['type'] = 's'
#              self.basis.append(bfn)
#            elif _type[g] == 'P':
#              bfn['type'] = 'px'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'py'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'pz'
#              self.basis.append(copy.deepcopy(bfn))
#            elif _type[g] == 'D':
#              bfn['type'] = 'dxx'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'dxy'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'dxz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'dyy'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'dyz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'dzz'
#              self.basis.append(copy.deepcopy(bfn))
#            elif _type[g] == 'F':
#              bfn['type'] = 'fxxx'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fxxy'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fxxz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fxyy'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fxyz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fxzz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fyyy'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fyyz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fyzz'
#              self.basis.append(copy.deepcopy(bfn))
#              bfn['type'] = 'fzzz'
#              self.basis.append(copy.deepcopy(bfn))
#            exp = []
#            cef = []
#          
#
#      movecs = stem+'.modat'
#      if os.path.exists(movecs):
#        self.getMO(movecs)
#
#  def getMO(self, mo_file):
#    """
#    setup self.n_ao
#          self.n_mo
#          self.occupation
#          self.mo_eigenvalues
#          self.mo_vectors
#          self.nuclear_repulsion
#    """
#    if os.path.exists(mo_file):
#      mo = open(mo_file, 'r')
#      mo_out = mo.readlines()
#      self.n_ao = int(mo_out[12].split( )[0])
#      self.n_mo = int(mo_out[13].split( )[0])
#      if self.n_ao % 3 != 0:
#        lines = self.n_ao / 3 + 1
#      else:
#        lines = self.n_ao / 3
#      self.occupation = \
#        [float(j) for i in mo_out[14:14+lines]\
#         for j in i.split( )]
#      self.mo_eigenvalues = \
#        [float(j) for i in mo_out[14+lines:14+(2*lines)]\
#         for j in i.split( )]
#
#      _mo = []
#      for k in range(self.n_mo):
#        start = 14+((2+k)*lines)
#        end = 14+((3+k)*lines)
#        vec = [float(j) for i in mo_out[start:end]\
#               for j in i.split( )]
#        _mo.append(vec)
#      self.mo_vectors = np.array(_mo)
#      self.nuclear_repulsion = float(mo_out[-1].split( )[1])
