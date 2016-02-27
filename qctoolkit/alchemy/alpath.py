import qctoolkit as qtk
import numpy as np
import copy
import universal as univ
from aljob import mutatePP

class AlPath(object):
  def __init__(self, ref, tar, **kwargs):
    self.setting = kwargs
    self.ref = univ.toInp(ref, **self.setting)
    self.tar = univ.toInp(tar, **self.setting)
    self.name = self.ref.molecule.name + '_' + self.tar.molecule.name
    self.mutation = {}

    # !!!!!!!!!!!!!!!!!!!!!!!
    # construct mutation list
    # tar_list: list of matching location
    ref = self.ref.molecule
    tar = self.tar.molecule
    tar_list = [True for i in range(tar.N)]
    to_void = True
    ref_ind = []
    tar_ind = []
    mutation_ind = [ref_ind, tar_ind]
    mutation_ref = []
    mutation_tar = []
    mutation_crd = []
    mutation_string = []
    # loop through every atom in reference system
    for i in range(ref.N):
      Zi = ref.type_list[i]
      Ri = ref.R[i]
      # loop through every atom in target system
      for j in range(tar.N):
        # every atom-pair distance
        Rj = tar.R[j]
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
        tar_ind.append(j)
        mutation_ref.append('VOID')
        mutation_tar.append(tar.type_list[j])
        mutation_crd.append(tar.R[j])
    for k in range(len(mutation_ref)):
      string = mutation_ref[k] + '2' + mutation_tar[k]
      mutation_string.append(string)
    self.mutation ={
                     'ind': mutation_ind,
                     'ref': mutation_ref,
                     'tar': mutation_tar,
                     'crd': mutation_crd,
                     'string': mutation_string,
                   }

  def __repr__(self):
    return self.name

  def write(self, name=None, **kwargs):
    if 'l' not in kwargs:
      kwargs['l'] = 0.2
    lambda_string = '_%03d' % (kwargs['l']*100)
    mol = self.ref.molecule.copy()
    mol.name = self.name + lambda_string
    N = mol.N
    mutation_list = copy.deepcopy(self.mutation['ind'][0])
    for i in self.mutation['ind'][1]:
      ind = N + i
      print ind, mol.N, i
      print self.mutation['tar']
      atom = self.mutation['tar'][ind]
      crd = self.mutation['crd'][ind]
      mol.addAtoms(atom, crd)
      mutation_list.append(ind)
    for string in set(self.mutation['string']):
      str_ind = [i for i in range(len(self.mutation['string']))\
                 if self.mutation['string'][i] == string]
      string = string + lambda_string
      mol.setAtoms(str_ind, string=string)

    qm_setting = {}
    if 'qm_setting' in self.setting:
      qm_setting.update(kwargs['qm_setting'])
    if 'qm_setting' in kwargs:
      qm_setting.update(kwargs['qm_setting'])
    if 'program' not in qm_setting:
      qm_setting['program'] = self.ref.setting['program']
    inp = qtk.QMInp(mol, **qm_setting)
    inp.write(name)
