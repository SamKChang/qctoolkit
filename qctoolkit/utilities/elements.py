import qctoolkit.data.elements as qel
import qctoolkit as qtk
import re

#################################
# element information utilities #
#################################
# load element data file one and for all
ve_list = qel.Elements.ve_list()
z_list = qel.Elements.z_list()
type_list = qel.Elements.type_list()
mass_list = qel.Elements.mass_list()

def n2ve(Zn):
  ref = re.sub('2[a-zA-Z].*','',Zn)
  tar = re.sub('.*[a-zA-Z]2','',Zn)
  tar = re.sub('_.*','',tar)
  # WARNING! symbol V is used for 
  if ref == 'V':
    ref = 'VOID'
    warning("VOID IS USED, symbol V is used for void "+\
            "instead of vanadium")
  if ve_list.has_key(Zn):
    return ve_list[Zn]
  elif ve_list.has_key(ref):
    return ve_list[ref]
  elif ve_list.has_key(tar):
    return ve_list[tar]
  else:
    qtk.exit("n2ve: element type " + Zn + " is not defined")

def Z2n(Z):
  if type_list.has_key(Z):
    return type_list[Z]
  else:
    qtk.exit("Z2n: atomic number " + str(Z) + " is not defined")
    #return Z
  
def n2Z(Zn):
  if z_list.has_key(Zn):
    return z_list[Zn]
  else:
    qtk.exit("n2Z: element type " + str(Zn) + " is not defined")

def n2Z0(Zn):
  if z_list.has_key(Zn):
    return z_list[Zn]
  else:
    return 0
  
def n2m(Zn):
  if mass_list.has_key(Zn):
    return mass_list[Zn]
  else:
    qtk.exit("n2Z: element type " + str(Zn) + " is not defined")

def qAtomName(query):
  if type(query) == str:
    if z_list.has_key(query):
      return str(query)
  elif type(query) == int or type(query) == float:
    if type_list.has_key(int(query)):
      return str(Z2n(query))
  else:
    qtk.exit("qAtom: element " + str(Zn) + " is not defined")

def qAtomicNumber(query):
  if type(query) == str:
    if z_list.has_key(query):
      return n2Z(query)
  elif type(query) == int or type(query) == float:
    if type_list.has_key(int(query)):
      return query
  else:
    qtk.exit("qAtom: element " + str(Zn) + " is not defined")

def isAtom(query):
  if z_list.has_key(query):
    return True
