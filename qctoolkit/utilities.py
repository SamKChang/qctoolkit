from math import pi ,sin, cos
import geometry as qg
import openbabel as ob
import numpy as np
import gc

def have_bond(obmol, type_a, type_b):
  bond = ob.OBBond()
  for i in range(obmol.NumBonds()):
    print type_a + "-" + type_b
    bond = obmol.GetBond(i)

  return True
  del bond

def qt2ob(qtmol):
  mol = ob.OBMol()
  for atom in xrange(qtmol.N):
    new_atom = mol.NewAtom()
    new_atom.SetAtomicNum(qtmol.Z[atom])
    new_atom.SetVector(qtmol.R[atom][0], 
                       qtmol.R[atom][1], 
                       qtmol.R[atom][2])
  mol.ConnectTheDots()
  return mol
  #del qtmol

def ob2qt(obmol):
  mol = qg.Molecule()
  mol.N = obmol.NumAtoms()
  atom = ob.OBAtom()
  bond = ob.OBBond()

  _Z = []
  _coordX = []
  _coordY = []
  _coordZ = []

  for i in range(1, obmol.NumAtoms()+1):
    atom = obmol.GetAtom(i)
    _Z.append(atom.GetAtomicNum())
    _coordX.append(atom.GetVector().GetX())
    _coordY.append(atom.GetVector().GetY())
    _coordZ.append(atom.GetVector().GetZ())
  for i in range(obmol.NumBonds()):
    bond = obmol.GetBond(0)
    #if bond is not None:
    obmol.DeleteBond(bond)
  for i in range(1, obmol.NumAtoms()+1):
    atom = obmol.GetAtom(1)
    obmol.DeleteAtom(atom)
  del obmol

  mol.Z = np.array(_Z)
  mol.R = np.vstack([_coordX, _coordY, _coordZ]).T

  return mol

  del _coordX, _coordY, coordZ, _Z, mol

def R(theta, u):
  return np.array(
    [[cos(theta) + u[0]**2 * (1-cos(theta)), 
      u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta), 
      u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
     [u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta),
      cos(theta) + u[1]**2 * (1-cos(theta)),
      u[1] * u[2] * (1 - cos(theta)) + u[0] * sin(theta)],
     [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
      u[1] * u[2] * (1-cos(theta)) - u[0] * sin(theta),
      cos(theta) + u[2]**2 * (1-cos(theta))]]
  )

def n2ve(Zn):
  ve_list = {
    'H' :  1,
    'He':  2,
    'Li':  3,
    'Be':  4,
    'B' :  3,
    'C' :  4,
    'N' :  5,
    'O' :  6,
    'F' :  7,
    'Ne':  8,
    'Na':  9,
    'Mg':  2,
    'Al':  3,
    'Si':  4,
    'P' :  5,
    'S' :  6,
    'Cl':  7,
    'Ar':  8,
    'K' :  9,
    'Ca':  2,
    'Ga':  3,
    'Ge':  4,
    'As':  5,
    'Se':  6,
  }
  if ve_list.has_key(Zn):
    return ve_list[Zn]
  else:
    return 0

def Z2n(Z):
  z_list = {
     1:'H' , 
     2:'He', 
     3:'Li', 
     4:'Be', 
     5:'B' , 
     6:'C' , 
     7:'N' , 
     8:'O' , 
     9:'F' , 
    10:'Ne', 
    11:'Na', 
    12:'Mg', 
    13:'Al', 
    14:'Si', 
    15:'P' , 
    16:'S' , 
    17:'Cl', 
    18:'Ar', 
    19:'K' , 
    20:'Ca', 
    31:'Ga', 
    32:'Ge', 
    33:'As', 
    34:'Se', 
  }
  if z_list.has_key(Z):
    return z_list[Z]
  else:
    return str(Z)
  
def n2Z(Zn):
  return {
    'H' :  1,
    'He':  2,
    'Li':  3,
    'Be':  4,
    'B' :  5,
    'C' :  6,
    'N' :  7,
    'O' :  8,
    'F' :  9,
    'Ne': 10,
    'Na': 11,
    'Mg': 12,
    'Al': 13,
    'Si': 14,
    'P' : 15,
    'S' : 16,
    'Cl': 17,
    'Ar': 18,
    'K' : 19,
    'Ca': 20,
    'Ga': 31,
    'Ge': 32,
    'As': 33,
    'Se': 34,
  }[Zn]
