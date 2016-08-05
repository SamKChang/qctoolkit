import numpy as np

def getPhi(inp, coords, **kwargs):
  if 'new' in kwargs and kwargs['new']:
    out = inp.mol.eval_gto("GTOval_sph", coords).T
  else:
    if inp._phi is None:
      inp._phi = inp.mol.eval_gto("GTOval_sph", coords).T
    out = inp._phi
  return out

def getDphi(inp, coords):
  if inp._dphi is None:
    inp._dphi = inp.mol.eval_gto("GTOval_ip_sph", coords, comp=3).T
  return inp._dphi

def getPsi(inp, coords=None, **kwargs):
  if coords is None:
    coords = inp.grid.points
  if 'dv' in kwargs:
    dv = kwargs['dv']
  else:
    dv = inp.dv
  if 'new' in kwargs and kwargs['new']:
    phi = getPhi(inp, coords, new=True)
    out = np.zeros(len(coords))
    for i in range(len(dv)):
      c_i = dv[i]
      phi_i = phi[i]
      out += c_i * phi_i
  else:
    if inp._psi is None:
      phi = getPhi(inp, coords)
      psi = np.zeros(len(coords))
      for i in range(len(dv)):
        c_i = dv[i]
        phi_i = phi[i]
        psi += c_i * phi_i
      inp._psi = psi
    out = inp._psi
  return out

def getDpsi(inp, coords=None, **kwargs):
  if inp._dpsi is None:
    if coords is None:
      coords = inp.grid.points
    if 'dv' in kwargs:
      dv = kwargs['dv']
    else:
      dv = inp.dv
    dphi = getDphi(inp, coords)
    dpsi = np.zeros([len(coords), 3])
    for i in range(len(dv)):
      c_i = dv[i]
      dphi_i = dphi[i]
      dpsi += c_i * dphi_i
    inp._dpsi = dpsi
  return inp._dpsi

def getRho(inp, coords=None, **kwargs):
  if coords is None:
    coords = inp.grid.points
  if 'new' in kwargs and kwargs['new']:
    out = getPsi(inp, coords, **kwargs)**2
  else:
    if inp._rho is None:
      inp._rho = getPsi(inp, coords, **kwargs)**2
    out = inp._rho
  return out

def getDrho(inp, coords=None, **kwargs):
  if inp._drho is None:
    if coords is None:
      coords = inp.grid.points
    psi = getPsi(inp, coords, **kwargs)
    dpsi = getDpsi(inp, coords, **kwargs)
    inp._drho = 2 * dpsi * psi[:, np.newaxis]
  return inp._drho

def getSigma(inp, coords=None, **kwargs):
  if inp._sigma is None:
    if coords is None:
      coords = inp.grid.points
    drho = getDrho(inp, coords, **kwargs)
    inp._sigma = np.sum(drho**2, axis=1)
  return inp._sigma

