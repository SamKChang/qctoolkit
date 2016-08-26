import grid_points as gp
import numpy as np
from libxc_dict import xc_dict
import libxc_interface as xcio

outer = np.outer
trace = np.trace
td = np.tensordot

def E(inp, coords = None, **kwargs):
  if coords is None:
    coords = inp.grid.points
  rho = gp.getRho(inp, coords)
  sigma = gp.getSigma(inp, coords)
  dft = inp.setting['dft_setting']
  if 'term' not in kwargs:
    e_k = np.zeros(len(coords))
    for fraction, kf in dft['K'].iteritems():
      e_k += fraction * xcio.exc(inp, kf, [rho, sigma], False)
    K = inp.grid.integrate(rho*e_k)
    V = trace(inp.dm.dot(inp.ext))
    int_vee_rho = td(inp.dm, inp.vee, axes=([0,1], [0,1]))
    U = trace(inp.dm.dot(int_vee_rho))/2.
    intc = xcio.exc(inp, dft['C'], [rho, sigma], False)
    intx = xcio.exc(inp, dft['X'], [rho, sigma], False)
    XC = inp.grid.integrate(rho * (intc + intx))
    E = K + V + U + XC
  elif kwargs['term'] == 'K':
    e_k = np.zeros(len(coords))
    for fraction, kf in dft['K'].iteritems():
      e_k += fraction * xcio.exc(inp, kf, [rho, sigma], False)
    E = inp.grid.integrate(rho*e_k)
  elif kwargs['term'] == 'vw':
    E = trace(inp.dm.dot(inp.kin))
  elif kwargs['term'] == 'U':
    int_vee_rho = td(inp.dm, inp.vee, axes=([0,1], [0,1]))
    E = trace(inp.dm.dot(int_vee_rho))/2.
  elif kwargs['term'] == 'V':
    E = trace(inp.dm.dot(inp.ext))
  elif kwargs['term'] == 'X':
    intx = xcio.exc(inp, dft['X'], [rho, sigma], False)
    E = inp.grid.integrate(rho * intx)
  elif kwargs['term'] == 'C':
    intc = xcio.exc(inp, dft['C'], [rho, sigma], False)
    E = inp.grid.integrate(rho * intc)
  elif kwargs['term'] in xc_dict.values() \
  or kwargs['term'] in xc_dict:
    epsilon = xcio.exc(inp, kwargs['term'], [rho, sigma])
    E = inp.grid.integrate(rho*epsilon)
  return E

def E_aufbau(inp, **kwargs):
  pass

def dE_ddv(inp, coords = None):
  if coords is None:
    coords = inp.grid.points

  # derivative is buggy... 

  rho = gp.getRho(inp, coords)
  sigma = gp.getSigma(inp, coords)
  dft = inp.setting['dft_setting']

  dEk_drho = np.zeros(len(coords))
  dEk_dsigma = np.zeros(len(coords))
  Ek = np.zeros(len(coords))
  for fraction, kf in dft['K'].iteritems():
    dEk_drho_f, dEk_dsigma_f = xcio.vxc(inp, kf, [rho, sigma], False)
    dEk_drho += fraction * dEk_drho_f
    dEk_dsigma += fraction * dEk_dsigma_f
    Ek += xcio.exc(inp, kf, [rho, sigma], False)

  dEc_drho, dEc_dsigma = xcio.vxc(inp, dft['C'], [rho, sigma], False)
  dEx_drho, dEx_dsigma = xcio.vxc(inp, dft['X'], [rho, sigma], False)
  #Ec = xcio.exc(inp, dft['C'], [rho, sigma], False)
  #Ex = xcio.exc(inp, dft['X'], [rho, sigma], False)
  dE_drho = dEk_drho + dEc_drho + dEx_drho
  dE_dsigma = dEk_dsigma + dEc_dsigma + dEx_dsigma
  #epsilon = Ek + Ec + Ex

  dE_kxc = np.zeros(len(inp.dv))
  for i in range(len(inp.dv)):
    drho_dci, dsigma_dci = drhoSigma_dc(inp, i)
    rho_dE_dc = dE_drho * drho_dci + dE_dsigma * dsigma_dci
    #E_drho_dc = drho_dci * epsilon
    integrand = rho_dE_dc
    dE_kxc[i] = inp.grid.integrate(integrand)

  vee_rho = td(inp.dm, inp.vee, axes=([0,1], [0,1]))
  dE_ee = 2 * vee_rho.dot(inp.dv)
  dE_ext = 2 * inp.ext.dot(inp.dv)
  out = dE_kxc + dE_ee + dE_ext

  return out

def drhoSigma_dc(inp, i, coords = None, **kwargs):
  if coords is None:
    coords = inp.grid.points

  phi = gp.getPhi(inp, coords, **kwargs)
  dphi = gp.getDphi(inp, coords)
  psi = gp.getPsi(inp, coords, **kwargs)
  dpsi = gp.getDpsi(inp, coords, **kwargs)
  drho = gp.getDrho(inp, coords, **kwargs)

  def drho_dc(i):
    return 2 * phi[i] * psi

  def nabla_drho_dc(i):
    term1 = 2 * dphi[i] * psi[:, np.newaxis]
    term2 = 2 * dpsi * phi[i][:, np.newaxis]
    return term1 + term2

  def dsigma_dc(i):
    return np.sum(2 * drho * nabla_drho_dc(i), axis=1)

  return drho_dc(i), dsigma_dc(i)

def eval_E(inp, dv):
  dv = inp.normalize(dv)
  inp.update(dv)
  coords = inp.grid.points
  dm = outer(dv, dv)
  rho = gp.getRho(inp, coords, dv=dv)
  sigma = gp.getSigma(inp, coords, dv=dv)
  dft = inp.setting['dft_setting']
  e_k = np.zeros(len(coords))
  for fraction, kf in dft['K'].iteritems():
    e_k += fraction * xcio.exc(inp, kf, [rho, sigma], False)
  K = inp.grid.integrate(rho*e_k)
  V = trace(dm.dot(inp.ext))
  int_vee_rho = td(dm, inp.vee, axes=([0,1], [0,1]))
  U = trace(dm.dot(int_vee_rho))/2.
  intc = xcio.exc(inp, dft['C'], [rho, sigma], False)
  intx = xcio.exc(inp, dft['X'], [rho, sigma], False)
  XC = inp.grid.integrate(rho * (intc + intx))
  E = K + V + U + XC
  inp.update(dv)
  return E

def E_dv(inp, dv):
  dv = inp.normalize(dv)
  coords = inp.grid.points
  dm = outer(dv, dv)
  rho = gp.getRho(inp, coords, dv=dv, new=True)
  sigma = gp.getSigma(inp, coords, dv=dv, new=True)
  dft = inp.setting['dft_setting']
  e_k = np.zeros(len(coords))
  for fraction, kf in dft['K'].iteritems():
    e_k += fraction * xcio.exc(inp, kf, [rho, sigma], False)
  K = inp.grid.integrate(rho*e_k)
  V = trace(dm.dot(inp.ext))
  int_vee_rho = td(dm, inp.vee, axes=([0,1], [0,1]))
  U = trace(dm.dot(int_vee_rho))/2.
  intc = xcio.exc(inp, dft['C'], [rho, sigma], False)
  intx = xcio.exc(inp, dft['X'], [rho, sigma], False)
  XC = inp.grid.integrate(rho * (intc + intx))
  E = K + V + U + XC
  return E

def eval_dE_ddv(inp, dv):
  coords = inp.grid.points
  dm = outer(dv,dv)

  rho = gp.getRho(inp, coords, dv=dv)
  sigma = gp.getSigma(inp, coords, dv=dv)
  dft = inp.setting['dft_setting']

  dEk_drho = np.zeros(len(coords))
  dEk_dsigma = np.zeros(len(coords))
  #Ek = np.zeros(len(coords))
  for fraction, kf in dft['K'].iteritems():
    dEk_drho_f, dEk_dsigma_f = xcio.vxc(inp, kf, [rho, sigma], False)
    dEk_drho += fraction * dEk_drho_f
    dEk_dsigma += fraction * dEk_dsigma_f
    #Ek += xcio.exc(inp, kf, [rho, sigma], False)

  dEc_drho, dEc_dsigma = xcio.vxc(inp, dft['C'], [rho, sigma], False)
  dEx_drho, dEx_dsigma = xcio.vxc(inp, dft['X'], [rho, sigma], False)
  #Ec = xcio.exc(inp, dft['C'], [rho, sigma], False)
  #Ex = xcio.exc(inp, dft['X'], [rho, sigma], False)
  dE_drho = dEk_drho + dEc_drho + dEx_drho
  dE_dsigma = dEk_dsigma + dEc_dsigma + dEx_dsigma
  #epsilon = Ek + Ec + Ex

  dE_kxc = np.zeros(len(dv))
  for i in range(len(dv)):
    drho_dci, dsigma_dci = drhoSigma_dc(inp, i, dv=dv)
    rho_dE_dc = dE_drho * drho_dci + dE_dsigma * dsigma_dci
    #E_drho_dc = drho_dci * epsilon
    integrand = rho_dE_dc# + E_drho_dc
    dE_kxc[i] = inp.grid.integrate(integrand)

  vee_rho = td(dm, inp.vee, axes=([0,1], [0,1]))
  dE_ee = 2 * vee_rho.dot(dv)
  dE_ext = 2 * inp.ext.dot(dv)
  out = dE_kxc + dE_ee + dE_ext

  return out
