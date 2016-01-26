#!/usr/bin/python

import qctoolkit as qtk
from qctoolkit.QM.atomicbasis_io import veMatrix
import glob
import numpy as np
import re
from scipy.interpolate import UnivariateSpline as uspl
from scipy.integrate import simps
import matplotlib.pyplot as plt

qtk.setting.quiet = True

# perturbation potentail
coord_i = np.array([
  [-0.966946,    0.549416,    0.000000],
  [-2.136715,    0.445692,    0.000000],
  [-0.335928,    1.725102,    0.000000],
  [0.000000,   -0.660843,    0.000000],
  [-0.405478,   -1.766520,    0.000000],
  [1.715445,   -0.252396,    0.000000]
])
coord_f = np.array([
  [-0.966946,    0.549416,    0.000000],
  [-2.136715,    0.445692,    0.000000],
  [-0.335928,    1.725102,    0.000000],
  [0.000000,   -0.660843,    0.000000],
  [-0.405478,   -1.766520,    0.000000],
  [1.715445,   -0.252396,    0.000000]
])
Z_i = [6, 8,9,6,8,17]
Z_f = [6, 8,17,6,8,9]

# setup variables
outs = sorted(glob.glob('*.log'))
Etot = []
Vn = []
dE = []
lmd = []

# load/calculate data
for o in outs:
  out = qtk.QMOut(o, program='gaussian')
  l = re.sub('OOFCl_0','',o)
  l = re.sub('.log','',l)
  lmd.append(float(l))
  Etot.append(out.Et)
  Vn.append(out.nuclear_repulsion)
  aovi = gcMatrix(out.basis, coord_i, Z_i)
  aovf = gcMatrix(out.basis, coord_f, Z_f)
  mov = np.dot(out.mo_vectors, np.dot(aovf-aovi, out.mo_vectors.T))
  print "diag.mov len ", len(np.diag(mov))
  print "occulation len ", len(out.occupation)
  dE_tmp = np.dot(np.diag(mov), out.occupation)
  dE.append(dE_tmp)
  """
  print mov
  print np.diag(mov)
print out.occupation
ind = [i for i in range(len(out.occupation)) if out.occupation[i] == 2.0]
print ind
print np.diag(mov)[ind[-1]+1]
"""
# convert to numpy array for post processing
Etot = np.array(Etot)
Vn = np.array(Vn)
dE = np.array(dE)
lmd = np.array(lmd)

# nuclear repulsion corretion
print Vn
dV = Vn[-1] - Vn[0]
Vcrr = np.array([Vn[0] + dV*i for i in lmd])
Ecrr = Etot - Vn + Vcrr
Ecrr_inter = uspl(lmd, Ecrr, k=4, s=0)
dEcrr = Ecrr_inter.derivative()
dE = dE + dV
dE_inter = uspl(lmd, dE, k=4, s=0)

print dEcrr.integral(lmd[0], lmd[-1])*627.509,
print " kcal/mol from finite difference"
print dE_inter.integral(lmd[0], lmd[-1])*627.509,
print " kcal/mol from analytic"
print (Etot[-1] - Etot[0])*627.509,
print " kcal/mol from nwchem"

# plot result
plt.figure('E_tot')
pCrr, = plt.plot(lmd, Ecrr_inter(lmd), label='Vn corrected')
pRow, = plt.plot(lmd, Etot, label='row data')
#plt.legend(handles=[pRow, pCrr])
plt.figure('dE_analytic')
anE, = plt.plot(lmd, dE, marker='s', color='b', markersize=8, label='analytic')
fdE, = plt.plot(lmd, dEcrr(lmd), marker='o', linestyle='--', color='r', label='finite difference')
#plt.legend(handles=[anE, fdE])
#plt.show()
