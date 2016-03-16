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
    qtk.warning("VOID IS USED, symbol V is used for void "+\
                "instead of vanadium")
  match = [m for m in ve_list.iterkeys() if m in Zn]
  match_tar = [m for m in ve_list.iterkeys() if m in tar]
  match_ref = [m for m in ve_list.iterkeys() if m in ref]
  if len(match) == 1:
    return ve_list[match[0]]
  elif len(match_ref) == 1:
    return ve_list[match_ref[0]]
  elif len(match_tar) == 1:
    return ve_list[match_tar[0]]
  else:
    qtk.exit("n2ve: element type " + Zn + " is not defined")

def Z2n(Z):
  if type_list.has_key(Z):
    return type_list[Z]
  elif type(Z) is float:
    return 'ATOM_%4.2f' % Z
    qtk.warning('Z2n: atomic number not defined, return HETATM')
  elif type(Z) is int:
    return 'ATOM_%d' % Z
    qtk.warning('Z2n: atomic number not defined, return HETATM')
  else:
    qtk.exit("Z2n: atomic number " + str(Z) + " is not defined")
    #return Z
  
def n2Z(Zn):
  match = [m for m in z_list.iterkeys() if m in Zn]
  if len(match) == 1:
    return float(z_list[match[0]])
  else:
    print Zn, 
    qtk.warning("n2Z: element type " + str(Zn) +\
                " is not defined, returning nuclear charge 0")
    return 0

def n2Z0(Zn):
  match = [m for m in z_list.iterkeys() if m in Zn]
  if len(match) == 1:
    return float(z_list[match[0]])
  else:
    return 0
  
def n2m(Zn):
  match = [m for m in mass_list.iterkeys() if m in Zn]
  if len(match) == 1:
    return float(mass_list[match[0]])
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
