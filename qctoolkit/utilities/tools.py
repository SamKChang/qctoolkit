from math import pi ,sin, cos
from qctoolkit.molecule import Molecule
import numpy as np
import qctoolkit as qtk
import yaml

def R(theta, u):
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

  sa = np.sqrt(1-ca**2) # sin(alpha) = sqrt(1-cos(alpha)^2)
  sc = np.sqrt(1-cc**2) # sin(gamma) = sqrt(1-cos(gamma)^2)
  cw = (cb - ca*cc)/(sa*sc)
  sw = np.sqrt(1-cw**2)

  return np.array(
    [
      [a * sc * sw, 0, 0],
      [a * cc, b, c * ca],
      [a * sc * cw, 0, c * sa],
    ]
  )

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
    'Eh': 1,
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
