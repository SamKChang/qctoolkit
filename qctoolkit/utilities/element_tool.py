import qctoolkit.data.elements as qel
import qctoolkit as qtk
import re
import numpy as np

#################################
# element information utilities #
#################################
# load element data file one and for all
ve_list = qel.Elements.ve_list()
z_list = qel.Elements.z_list()
type_list = qel.Elements.type_list()
mass_list = qel.Elements.mass_list()

def e2ve(e):
  pattern = re.compile('\d[spdf]')
  ve = 0
  ecfg_lst = e.eleconfig.split()
  print ecfg_lst
  ecfg_str = filter(pattern.match, ecfg_lst)
  for s in ecfg_str[-2:]:
    n = re.sub(pattern, '', s)
    try:
      ve = ve + int(n)
    except ValueError:
      ve = ve + 1

  if e.period > 2:
    if e.group == 18:
      ve = 8
    elif e.group < 3:
      ve = ve + 8
    else:
      if e.period == 4:
        if e.group < 11:
          ve = ve + 8
      elif e.period == 5:
        if e.group < 8:
          ve = ve + 8
      elif e.period == 6:
        if e.group < 5:
          ve = ve + 8
    
  elif e.period > 1:
    if e.group < 13:
      ve = ve + 2

  return ve

def n2ve(Zn):
  ref = re.sub('2[a-zA-Z].*','',Zn)
  tar = re.sub('.*[a-zA-Z]2','',Zn)
  tar = re.sub('_.*','',tar)
  # WARNING! symbol V is used for 
  match = [m for m in ve_list.iterkeys() if m in Zn]
  match_tar = [m for m in ve_list.iterkeys() if m in tar]
  match_ref = [m for m in ve_list.iterkeys() if m in ref]
  mlen = [len(s) for s in match]
  rlen = [len(s) for s in match_ref]
  tlen = [len(s) for s in match_tar]
  if len(match) > 0:
    ind = np.argmax(mlen)
    return ve_list[match[ind]]
  elif len(match_ref) > 0:
    ind = np.argmax(rlen)
    return ve_list[match_ref[ind]]
  elif len(match_tar) > 0:
    ind = np.argmax(tlen)
    return ve_list[match_tar[ind]]
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
  mlen = [len(s) for s in match]
  if len(match) > 0:
    ind = np.argmax(mlen)
    return float(z_list[match[ind]])
  else:
    qtk.warning("n2Z: element type " + str(Zn) +\
                " is not defined, returning nuclear charge 0")
    return 0

def n2Z0(Zn):
  match = [m for m in z_list.iterkeys() if m in Zn]
  mlen = [len(s) for s in match]
  if len(match) > 0:
    ind = np.argmax(mlen)
    return float(z_list[match[ind]])
  else:
    return 0
  
def n2m(Zn):
  match = [m for m in mass_list.iterkeys() if m in Zn]
  mlen = [len(s) for s in match]
  if len(match) > 0:
    ind = np.argmax(mlen)
    return float(mass_list[match[ind]])
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
