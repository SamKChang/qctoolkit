import qctoolkit as qtk
import numpy as np
import copy
import universal as univ
from aljob import mutatePP

class AlPath(object):
  """
  alchemical path object to interpolate pseudo potentail and 
  perform lambda scan of E(lambda)
  """
  def __init__(self, ref, tar, **kwargs):
    self.setting = kwargs
    if 'pp_string' not in kwargs:
      self.setting['pp_string'] = 'pbe'
    self.ref = univ.toInp(ref, **self.setting)
    self.tar = univ.toInp(tar, **self.setting)
    self.name = self.ref.molecule.name + '_' + self.tar.molecule.name
    self.mutation = {}

    qm_setting = {}
    if 'qm_setting' in self.setting:
      qm_setting.update(kwargs['qm_setting'])
    if 'qm_setting' in kwargs:
      qm_setting.update(kwargs['qm_setting'])
    if 'program' not in qm_setting:
      qm_setting['program'] = self.ref.setting['program']
    self.qm_setting = qm_setting

    # !!!!!!!!!!!!!!!!!!!!!!!
    # construct mutation list
    # tar_list: list of matching location
    ref = self.ref.molecule
    tar = self.tar.molecule
    tar_list = [True for i in range(tar.N)]
    to_void = True
    ref_ind = []
    tar_ind = []
    mutation_fix = []
    mutation_ind = ref_ind
    mutation_ref = []
    mutation_tar = []
    mutation_crd = []
    mutation_string = []
    # loop through every atom in reference system
    for i in range(ref.N):
      Zi = ref.type_list[i]
      Ri = copy.deepcopy(ref.R[i])
      # loop through every atom in target system
      for j in range(tar.N):
        # every atom-pair distance
        Rj = copy.deepcopy(tar.R[j])
        diff = np.linalg.norm(Ri-Rj)
        # search for atoms at same location
        if diff < 10E-6:
          to_void = False
          tar_list[j] = False
          Zj = tar.type_list[j]
          if Zi != Zj:
            ref_ind.append(i)
            mutation_ref.append(Zi)
            mutation_tar.append(Zj)
            mutation_crd.append(Ri)
          else:
            mutation_fix.append(i)
      # when no distance matched, set atom target to void
      if to_void:
        ref_ind.append(i)
        mutation_ref.append(Zi)
        mutation_tar.append('VOID')
        mutation_crd.append(Ri)
    # loop through every target atom for non-matching atoms
    N = ref.N
    for j in range(tar.N):
      if tar_list[j]:
        ref_ind.append(j)
        mutation_ref.append('VOID')
        mutation_tar.append(tar.type_list[j])
        mutation_crd.append(tar.R[j])
    for k in range(len(mutation_ref)):
      string = mutation_ref[k] + '2' + mutation_tar[k]
      if self.setting['pp_string']:
        string = string + '_' + self.setting['pp_string']
      mutation_string.append(string)
    self.mutation ={
                     'ind': mutation_ind,
                     'ref': mutation_ref,
                     'tar': mutation_tar,
                     'crd': mutation_crd,
                     'fix': mutation_fix,
                     'string': mutation_string,
                   }


  def __repr__(self):
    return self.name

  def setLambda(self, **kwargs):
    if 'l' not in kwargs:
      kwargs['l'] = 0.2
    # dirty fix, grep '_[0-9]' pattern in QM.write function
    lambda_string = '_%03d' % (kwargs['l']*100)
    del kwargs['l']
    mol = self.ref.molecule.copy()
    mol.name = self.name + lambda_string
    N = mol.N
    string_dict = {}
    for i in range(len(self.mutation['ind'])):
      if i >= N:
        mol.addAtoms(self.mutation['tar'][i], 
                     self.mutation['crd'][i])
      ind = self.mutation['ind'][i]
      pp_string = self.mutation['string'][i] + lambda_string
      if pp_string not in string_dict:
        string_dict[pp_string] = [ind]
      else:
        string_dict[pp_string].append(ind)

    itr = 1
    for pp_string, ind in string_dict.iteritems():
      mol.setAtoms(ind, string=pp_string)
      Zn = mol.type_list[ind[0]]
      mol.setAtoms(ind, element = Zn + str(itr))
      itr = itr + 1

    return mol

  def write(self, name=None, **kwargs):
    molecule = self.setLambda(**kwargs)
    self.qm_setting.update(kwargs)
    inp = qtk.QMInp(molecule, **self.qm_setting)
    inp.write(name)

  def getInp(self, **kwargs):
    molecule = self.setLambda(**kwargs)
    self.qm_setting.update(kwargs)
    inp = qtk.QMInp(molecule, **self.qm_setting)
    return inp

  def writeAll(self, name, **kwargs):
    if 'dl' not in kwargs:
      dl = 0.2
    if 'l' in kwargs: 
      del kwargs['l']
    for l in np.arange(0,1,dl):
      l_str = '_%03d' % (l * 100)
      new_name = name + l_str
      kwargs['l'] = l
      self.write(new_name, **kwargs)

  def run(self, name=None, **kwargs):
    molecule = self.setLambda(**kwargs)
    self.qm_setting.update(kwargs)
    inp = qtk.QMInp(molecule, **self.qm_setting)
    return inp.run(name)

  def runAll(self, name=None, **kwargs):
    if 'dl' not in kwargs:
      kwargs['dl'] = 0.2
    if 'l' in kwargs: 
      del kwargs['l']
    out = []
    for l in np.arange(0,1.01,kwargs['dl']):
      l_str = '_%03d' % (l * 100)
      new_name = name + l_str
      kwargs['l'] = l
      out.append(self.run(new_name, **kwargs))
    return out
