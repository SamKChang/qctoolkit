import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import AtomicBasisInput
from qctoolkit.QM.atomicbasis_io import AtomicBasisOutput
import os, sys, copy, shutil, re
import numpy as np
import qctoolkit.QM.qmjob as qmjob

class inp(AtomicBasisInput):
  def __init__(self, molecule, **kwargs):
    AtomicBasisInput.__init__(self, molecule, **kwargs)
    self.setting.update(kwargs)

  def run(self, name=None, **kwargs):
    # creat_folder routine check default name
    name = self.create_folder(name)
    cwd = os.getcwd()
    os.chdir(name)
    inp = name + '.inp'
    self.write(inp)
    try:
      out = qmjob.QMRun(inp, 'nwchem', **kwargs)
    except:
      qtk.warning("qmjob finished unexpectedly for '" + \
                  name + "'")
      out = AtomicBasisOutput(program='nwchem')
    finally:
      os.chdir(cwd)
    return out
    
  def write(self, name=None):
    molecule = copy.deepcopy(self.molecule)
    self.cm_check(molecule)
    if name:
      nameStr = re.sub('\.inp', '', name)
      out = os.path.splitext(name)
      if out[1] != '.inp':
        name = out[0] + '.inp'
      if os.path.exists(name) and not qtk.setting.no_warning:
        qtk.prompt(name + ' exists, overwrite?')
        try:
          os.remove(name)
        except:
          qtk.exit("can not remove file: " + name)
    else:
      nameStr = re.sub('\.inp', '', molecule.name)
    inp = sys.stdout if not name else open(name, 'w')

    # mode setup
    if self.setting['mode'] == 'single_point':
      operation = 'energy'
    elif self.setting['mode'] == 'geopt':
      operation = 'optimize'

    # Theory setupt
    dft = {
      'pbe': 'xpbe96 cpbe96',
      'pbe0': 'pbe0',
      'blyp': 'becke88 lyp',
      'b3lyp': 'b3lyp',
      'bp91': 'becke88 perdew91',
      'bp86': 'becke88 perdew86',
      'pw91': 'xperdew91 perdew91',
    }
    others = {
      'rhf', 'rohf', 'uhf',
      'mp2', 'ccsd', 'ccsd(t)',
      'md'
    }
    if self.setting['theory'] in dft:
      module = 'dft'
      xc = dft[self.setting['theory']]
    elif self.setting['theory'] in others:
      module = self.setting['theory']

    # print file
    inp.write('title %s\n' % nameStr)
    inp.write('start %s\n' % nameStr)
    inp.write('echo\n')
    if self.setting['fix_molecule']\
    and self.setting['fix_molecule']:
      fix_mol = 'nocenter'
    else:
      fix_mol = ''
    inp.write('\ngeometry units angstrom %s\n' % fix_mol)
    if fix_mol:
      inp.write('symmetry group c1\n')
    if 'nuclear_charges' in self.setting:
      new_Z = self.setting['nuclear_charges']
    else:
      new_Z = []
    z_itr = 0
    for i in range(molecule.N):
      n_charge = ''
      e_str = molecule.type_list[i]
      if new_Z and z_itr < len(new_Z):
        if new_Z[z_itr][0] == i+1:
          n_charge = 'charge %.3f' % float(new_Z[z_itr][1])
          e_str = new_Z[z_itr][2]
          molecule.string[i] = e_str
          z_itr = z_itr + 1
      inp.write(' %-8s % 8.4f % 8.4f % 8.4f %s\n' % \
                (e_str,
                 molecule.R[i, 0],
                 molecule.R[i, 1],
                 molecule.R[i, 2],
                 n_charge
               ))
    inp.write('end\n\n')
    inp.write('basis\n')
    if self.setting['basis_set'] != 'gen':
      eStr = []
      for e in range(len(molecule.Z)):
        if molecule.string[e]:
          eString = molecule.string[e]
        else:
          eString = molecule.type_list[e]
        if eString not in set(eStr):
          eStr.append(eString)
          inp.write(' %-8s library %2s %s\n' % (\
            eString,
            molecule.type_list[e],
            self.setting['basis_set'].upper()),
          )
    inp.write('end\n\n')
    if module == 'dft':
      inp.write('dft\n')
      inp.write(' xc %s\n' % xc)
      inp.write('end\n\n')
    elif module == 'scf':
      pass
    if molecule.charge != 0:
      inp.write('charge % 5.2f\n' % molecule.charge)
    if molecule.multiplicity != 1:
      inp.write('multiplicity %d\n' % molecule.multiplicity)

    inp.write('\ntask %s %s\n' % (module, operation))

    if name:
      inp.close()


class out(AtomicBasisOutput):
  """
  directly parse vasp xml output, 'vasprun.xml'
  converged energy, system info, and scf steps are extracted
  """
  def __init__(self, qmout=None, **kwargs):
    AtomicBasisOutput.__init__(self, qmout=None, **kwargs)
    if qmout:
      stem, ext = os.path.splitext(qmout)
      outfile = open(qmout, 'r')
      data = outfile.readlines()
      Et = filter(lambda x: 'Total' in x and 'energy' in x, data)
      try:
        self.Et = float(Et[-1].split( )[-1])
      except:
        self.Et = np.nan
      n_basis = filter(lambda x: 'number of functions' in x, data)
      try:
        self.n_basis = int(n_basis[-1].split(':')[1])
      except:
        self.n_basis = np.nan

      ######################################
      # extract basis function information #
      ######################################
      basis_P = re.compile(r"  [0-9] [A-Z] +")
      batom_P = re.compile(r"^  [A-Z][a-z\ ] *\(")
      coord_P = re.compile(r"^ [A-Z][a-z ] +[-0-9\.]" + \
                           r"{10,}[-0-9\. ]+$")
      basisStr = filter(basis_P.match, data)
      batomStr = filter(batom_P.match, data)
      coordStr = filter(coord_P.match, data)

      _exponents = [float(filter(None, s.split(' '))\
        [2]) for s in basisStr]
      _coefficients = [float(filter(None, s.split(' '))\
        [3]) for s in basisStr]
      _N = [int(filter(None, s.split(' '))[0])\
        for s in basisStr]
      _type = [filter(None, s.split(' '))[1]\
        for s in basisStr]
      ao_keys = np.diff(_N)
      self.ao_keys = [i+1 for i in range(len(ao_keys))\
        if ao_keys[i] < 0]
      self.ao_keys.append(len(_N))
      self.ao_keys.insert(0, 0)
      _atoms = [filter(None, s.split(' '))[0]\
        for s in batomStr]
      self.type_list = [filter(None, s.split(' '))[0]\
        for s in coordStr]
      self.R = np.array([filter(None, s.split(' '))[1:4]\
        for s in coordStr]).astype(float)
      self.R_bohr = 1.889725989 * self.R

      _N.append(0)
      self.basis = []
      print _N
      print len(_N)
      print self.ao_keys
      for i in range(len(self.type_list)):
        e = self.type_list[i]
        center = self.R_bohr[i]
        ind = _atoms.index(e)
        bfn_base = {}
        bfn_base['atom'] = e
        bfn_base['center'] = center
        bfn_base['index'] = i
        exp = []
        cef = []
        for g in range(self.ao_keys[ind], self.ao_keys[ind+1]):
          exp.append(_exponents[g])
          cef.append(_coefficients[g])
          print _N[g], g
          if _N[g] != _N[g+1]:
            bfn = copy.deepcopy(bfn_base)
            bfn['exponents'] = copy.deepcopy(exp)
            bfn['coefficients'] = copy.deepcopy(cef)
            if _type[g] == 'S':
              bfn['type'] = 's'
              self.basis.append(bfn)
            elif _type[g] == 'P':
              bfn['type'] = 'px'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'py'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'pz'
              self.basis.append(copy.deepcopy(bfn))
            elif _type[g] == 'D':
              bfn['type'] = 'dxx'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'dxy'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'dxz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'dyy'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'dyz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'dzz'
              self.basis.append(copy.deepcopy(bfn))
            elif _type[g] == 'F':
              bfn['type'] = 'fxxx'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fxxy'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fxxz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fxyy'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fxyz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fxzz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fyyy'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fyyz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fyzz'
              self.basis.append(copy.deepcopy(bfn))
              bfn['type'] = 'fzzz'
              self.basis.append(copy.deepcopy(bfn))
            exp = []
            cef = []
          

      movecs = stem+'.modat'
      self.getMO(movecs)

  def getMO(self, mo_file):
    """
    setup self.n_basis
          self.occupation
          self.eigenvalues
          self.mo
          self.nuclear_repulsion
    """
    if os.path.exists(mo_file):
      mo = open(mo_file, 'r')
      mo_out = mo.readlines()
      #self.n_basis = int(mo_out[12].split( )[0])
      self.n_mo = int(mo_out[13].split( )[0])
      lines = self.n_basis / 3 + 1
      self.occupation = \
        [float(j) for i in mo_out[14:14+lines]\
         for j in i.split( )]
      self.mo_eigenvalues = \
        [float(j) for i in mo_out[14+lines:14+(2*lines)]\
         for j in i.split( )]

      _mo = []
      for k in range(self.n_mo):
        start = 14+((2+k)*lines)
        end = 14+((3+k)*lines)
        vec = [float(j) for i in mo_out[start:end]\
               for j in i.split( )]
        _mo.append(vec)
      self.mo = np.array(_mo)
      self.nuclear_repulsion = float(mo_out[-1].split( )[1])
