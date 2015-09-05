#!/usr/bin/python

import qctoolkit as qtk
import qctoolkit.alchemy as qal
import glob, copy, os, re, shutil
import numpy as np

if os.path.exists('inp'):
  shutil.rmtree('inp')
os.makedirs('inp')

dopt = [[ False, 2.1163, 2.2244], 
        [1.5190,  False, 1.9307],
        [1.5677, 1.8524,  False]]

xyz_file = glob.glob('structure/*.xyz')
inp = qtk.QMInp(xyz_file[0],'cpmd')
inp.setCelldm([20,15,15,0,0,0])
inp.setShift([7.5,7.5,7.5])
inp.saveDensity()
ref_list = []
tar_list = []
name_list = []
for xyz in xyz_file:
  mol1 = qtk.Molecule()
  mol1.read(xyz)
  ref_list.append(mol1)

  atom, coord = mol1.extract(1)
  tars = []
  for d_tar in [1, 1.5, 2, 2.5]:
    mol2 = qtk.Molecule()
    mol2.read(xyz)
    mol2.stretch(1, [2,1], d_tar-coord[0])
    tars.append(mol2)
  tar_list.append(tars)

  stem, ext = os.path.splitext(xyz)
  stem = re.sub('.*\/','',stem)
  name_list.append(re.sub('-eq', '', stem))

def verticalMutation(ref, tar, dopt):
  d_r = ref.extract(1)[1][0]
  d_t = tar.extract(1)[1][0]
  ref_ver = copy.deepcopy(tar)
  ref_ver.stretch(1,[2,1],d_r - d_t)
  path_ref = qal.PathScan(ref, ref_ver, inp_base=inp)
  tar_ver = copy.deepcopy(ref)
  tar_ver.stretch(1,[2,1],d_t - d_r)
  path_tar = qal.PathScan(tar_ver, tar, inp_base=inp)
  ref_opt = copy.deepcopy(ref)
  ref_opt.stretch(1,[2,1],dopt - d_r)
  tar_opt = copy.deepcopy(tar)
  tar_opt.stretch(1,[2,1],dopt - d_t)
  path_opt = qal.PathScan(ref_opt, tar_opt, inp_base=inp)
  return path_ref, path_tar, path_opt

for i in range(len(ref_list)):
  for f in range(len(tar_list)):
    if i != f:
      mol_i = ref_list[i]
      d_o = dopt[i][f]
      itr = 0
      for mol_f in tar_list[f]:
        path_ref, path_tar, path_opt = \
          verticalMutation(mol_i, mol_f, d_o)
        d_name = "_d%02d" % int(mol_f.extract(1)[1][0]*10)
        for l in np.linspace(0,1,21):
          l_str = "_l%03d" % (l*100)
          name = 'inp/' + name_list[i] + '-' + name_list[f]
          rname = name+"_ref"+l_str
          tname = name+d_name+l_str
          oname = name+"_opt"+l_str
          if itr==0:
            inp = path_ref.l(l)
            inp.setInfo(rname)
            inp.write(rname+'.inp')

            inp = path_opt.l(l)
            inp.setInfo(oname)
            inp.write(oname+'.inp')

  
          inp = path_tar.l(l)
          inp.setInfo(tname)
          inp.write(tname+'.inp')
        itr += 1
  
