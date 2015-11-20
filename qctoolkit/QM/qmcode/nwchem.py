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
      for e in range(len(molecule.Z)):
        if molecule.string[e]:
          eString = molecule.string[e]
        else:
          eString = molecule.type_list[e]
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
      nbasis = filter(lambda x: 'number of functions' in x, data)
      try:
        self.nbasis = int(nbasis[-1].split(':')[1])
      except:
        self.nbasis = np.nan

  def getMO(self, mo_file):
    """
    return self.n_basis
           self.occupation
           self.eigenvalues
           self.mo
           self.nuclear_repulsion
    """
    if os.path.exists(mo_file):
      mo = open(mo_file, 'r')
      mo_out = mo.readlines()
      self.n_basis = int(mo_out[13].split( )[0])
      lines = self.n_basis / 3 + 1
      self.occupation = \
        [float(j) for i in mo_out[14:14+lines]\
         for j in i.split( )]
      self.mo_eigenvalues = \
        [float(j) for i in mo_out[14+lines:14+(2*lines)]\
         for j in i.split( )]

      _mo = []
      for k in range(self.n_basis):
        start = 14+((2+k)*lines)
        end = 14+((3+k)*lines)
        vec = [float(j) for i in mo_out[start:end]\
               for j in i.split( )]
        _mo.append(vec)
      self.mo = np.array(_mo)
      self.nuclear_repulsion = float(mo_out[-1].split( )[1])
