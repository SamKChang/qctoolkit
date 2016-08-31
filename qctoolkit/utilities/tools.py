from math import pi ,sin, cos, ceil
from qctoolkit.molecule import Molecule
import numpy as np
import qctoolkit as qtk
import yaml
import pickle

def flatten(container):
  for i in container:
    if isinstance(i, (list,tuple)):
      for j in flatten(i):
        yield j
    else:
      yield i

def pdump(obj, name):
  with open(name, 'wb') as pfile:
    pickle.dump(obj, pfile)

def pload(name):
  with open(name, 'rb') as pfile:
    return pickle.load(pfile)

def CoulombMatrix(molecule, dim=None):
  mol = toMolecule(molecule)
  if dim is None:
    dim = mol.N
  if dim < mol.N:
    qtk.exit("Coulomb matrix dimension must greater than" +\
             " the number of atoms")
  M = np.zeros((dim, dim))
  for i in range(mol.N):
    for j in range(i, mol.N):
      if i == j:
        M[i, j] = 0.5 * mol.Z[i] ** 2.4
      else:
        Rij = np.linalg.norm(mol.R[i] - mol.R[j])
        M[i, j] = mol.Z[i] * mol.Z[j] / Rij
        M[j, i] = M[i, j]
  order = list(np.argsort(sum(M ** 2)))
  order.reverse()
  out = M[:, order]
  return out[order, :]

def R(theta, u):
  u = np.array(u) / np.linalg.norm(u)
  return np.array(
    [[cos(theta) + u[0]**2 * (1-cos(theta)), 
      u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta), 
      u[0] * u[2] * (1-cos(theta)) + u[1] * sin(theta)],
     [u[0] * u[1] * (1-cos(theta)) + u[2] * sin(theta),
      cos(theta) + u[1]**2 * (1-cos(theta)),
      u[1] * u[2] * (1-cos(theta)) - u[0] * sin(theta)],
     [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
      u[1] * u[2] * (1-cos(theta)) + u[0] * sin(theta),
      cos(theta) + u[2]**2 * (1-cos(theta))]]
  )

def fractionalMatrix(celldm_list):
  a = celldm_list[0]
  b = celldm_list[1]
  c = celldm_list[2]
  ca = celldm_list[3] # cos(alpha)
  cb = celldm_list[4] # cos(beta)
  cc = celldm_list[5] # cos(gamma)
  #sa = np.sqrt(1-ca**2) # sin(alpha) = sqrt(1-cos(alpha)^2)
  sc = np.sqrt(1-cc**2) # sin(gamma) = sqrt(1-cos(gamma)^2)
  #cw = (cb - ca*cc)/(sa*sc)
  #sw = np.sqrt(1-cw**2)
  v = np.sqrt(1 - ca**2 - cb**2 - cc**2 + 2*ca*cb*cc)

  # fractional transformation from wiki
  return np.array(
    [
      [a, b*cc, c*cb],
      [0, b*sc, c*(ca - cb*cc)/sc],
      [0, 0, c*v/sc],
    ]
  )

def cellVec2celldm(v):
  a = np.linalg.norm(v[0])
  b = np.linalg.norm(v[1])
  c = np.linalg.norm(v[2])
  cosA = np.dot(v[1], v[2]) / (b*c)
  cosB = np.dot(v[0], v[2]) / (a*c)
  cosC = np.dot(v[0], v[1]) / (a*b)
  return [a, b, c, cosA, cosB, cosC]

#  return np.array(
#    [
#      [a * sc * sw, 0, 0],
#      [a * cc, b, c * ca],
#      [a * sc * cw, 0, c * sa],
#    ]
#  )

def unscaledCelldm(celldm, scale):
  celldm_new = [celldm[i]/scale[i] for i in range(3)]
  angle = celldm[3:]
  celldm_new.extend(angle)
  return celldm_new
  
def lattice2celldm(lattice):
  celldm = [round(np.linalg.norm(a), 8) for a in lattice]
  for i in range(3):
    j = (i + 1) % 3
    k = (i + 2) % 3
    vj = lattice[j] / celldm[j]
    vk = lattice[k] / celldm[k]
    celldm.append(round(np.dot(vj, vk), 8))
  return celldm

def celldm2lattice(celldm, **kwargs):
  scale = [1.0 for i in range(3)]
  if 'scale' in kwargs:
    scale = kwargs['scale']
  celldm = unscaledCelldm(celldm, scale)
  fm = fractionalMatrix(celldm)
  return np.dot(fm, np.eye(3)).T

def fractional2xyz(R_scale, lattice_celldm):
  scale = [ceil(i) for i in np.max(R_scale,axis=0)]
  dim = np.array(lattice_celldm).shape
  if len(dim) == 2:
    assert dim == (3, 3)
    lattice = np.array(lattice_celldm) / np.array(scale)
  else:
    assert dim == (6,)
    lattice = celldm2lattice(lattice_celldm) / np.array(scale)
  R = np.dot(R_scale, lattice)
  return np.array(R)

#def fractional2xyz(R_scale, celldm, scale=[1,1,1]):
#  celldm = unscaledCelldm(celldm, scale)
#  fm = fractionalMatrix(celldm)
#  return np.dot(fm, R_scale.T).T

def xyz2fractional(R, celldm, scale=[1,1,1]):
  celldm = unscaledCelldm(celldm, scale)
  fm = fractionalMatrix(celldm)
  return np.dot(np.linalg.inv(fm), R.T).T

def convE(source, units, separator=None):
  def returnError(ioStr, unitStr):
    msg = 'supported units are:\n'
    for key in Eh.iterkeys():
      msg = msg + key + '\n'
    qtk.report(msg, color=None)
    qtk.exit(ioStr + " unit: " + unitStr + " is not reconized")

  EhKey = {
    'ha': 'Eh',
    'eh': 'Eh',
    'hartree': 'Eh',
    'ry': 'Ry',
    'j': 'J',
    'joule': 'J',
    'kj/mol': 'kJ/mol',
    'kjmol': 'kJ/mol',
    'kjm': 'kJ/mol',
    'kj': 'kJ/mol', # assume no kilo Joule!
    'kcal/mol': 'kcal/mol',
    'kcalmol': 'kcal/mol',
    'kcal': 'kcal/mol', # assume no kcal!
    'kcm': 'kcal/mol',
    'ev': 'eV',
    'cminv': 'cmInv',
    'cminverse': 'cmInv',
    'icm': 'cmInv',
    'cm-1': 'cmInv',
    'k': 'K',
    'kelvin': 'K',
  }
 
  Eh = {
    'Eh': 1.0,
    'Ry': 2.0,
    'eV': 27.211396132,
    'kcal/mol': 627.509469,
    'cmInv': 219474.6313705,
    'K': 3.15774646E5,
    'J': 4.3597443419E-18,
    'kJ/mol': 2625.49962
  }

  if not separator:
    separator='-'
  unit = units.split(separator)
  if len(unit) != 2:
    qtk.exit("problem with unit separator '%s'" % separator)
  if unit[0].lower() != 'hartree' and unit[0].lower() != 'eh':
    if unit[0].lower() in EhKey:
      unit0 = EhKey[unit[0].lower()]
      source = source / Eh[unit0]
    else: returnError('input', unit[0])
  if unit[1].lower() not in EhKey: 
    returnError('output', unit[1])
  else:
    unit1 = EhKey[unit[1].lower()]
  return source * Eh[unit1], unit1

def imported(module):
  try:
    __import__(module)
  except ImportError:
    return False
  else:
    return True

def toMolecule(input_data, **kwargs):
  if type(input_data) is not Molecule:
    try:
      return Molecule(input_data, **kwargs)
    except:
      pass
  else:
    return input_data

def numberToBase(n, b):
  if n == 0:
    return [0]
  digits = []
  while n:
    digits.append(int(n % b))
    n /= b
  return digits[::-1]

def partialSum(iterable):
  total = 0
  for i in iterable:
    total += i
    yield total

def listShape(input_list):
  if type(input_list) == list:
    if type(input_list[0]) != list:
      return len(input_list)
    else:
      return [listShape(sublist) for sublist in input_list]
