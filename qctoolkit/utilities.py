from math import pi ,sin, cos
import geometry as qg
#import openbabel as ob
import numpy as np
#import gc
import fileinput
import sys
import inspect
import multiprocessing as mp
import operator
from compiler.ast import flatten

def parallelize(target_function, 
                input_list, 
                threads,
                **kwargs):

  if 'block_size' in kwargs:
    block_size = kwargs['block_size']
  else:
    block_size = len(input_list)/(threads*3)

  #############################################
  # runing target function of a single thread #
  #############################################
  def run_jobs(q_in, q_out):
    for inps in iter(q_in.get, None):
      ind = inps[-1]    # index of job
      inps = inps[:-1]  # actual input sequence
      out = []
      for args in inps:
        if type(args[-1]) == dict: # check known args input
          kwargs = args[-1]
          args = args[:-1]  
          out.append(target_function(*args, **kwargs))
        else:
          out.append(target_function(*args))
      if out != None:
        q_out.put([out, ind]) # output result with index
      q_in.task_done()
    q_in.task_done() # task done for 'None' if q_in finished
  ###### end of single thread definition ######

  # devide input_list into chunks according to block_size
  def chunks(_list, _size):
    for i in xrange(0, len(_list), _size):
      yield _list[i:i+_size]
  input_block = list(chunks(input_list, block_size))

  # setup empty queue
  output_stack = []
  output = []
  qinp = mp.JoinableQueue()
  qout = mp.Queue()

  # start process with empty queue
  for thread in range(threads):
    p =  mp.Process(target=run_jobs, args=(qinp, qout))
    p.daemon = True # necessary for terminating finished thread
    p.start()

  # put I/O data into queue for parallel processing
  index = range(len(input_block))
  for ind, inps in zip(index, input_block):
    inps.append(ind) # append inp index
    qinp.put(inps)   # put inp to input queue
  qinp.join()       # wait for jobs to finish

  # 'while not queue.empty' is NOT reliable
  if not qout.empty():
    for i in range(len(input_block)):
      output_stack.append(qout.get())

  if len(output_stack)>0:
    # sort/restructure output according to input order
    output_stack = sorted(output_stack, key=operator.itemgetter(1))
    # loop though all input for corresponding output
    for data_out in output_stack: 
      # if output is list of class, in-line iteration doesn't work
      output.append(data_out[0])
    return flatten(output)

class bcolors:
  HEADER = '\033[95m'
  OKBLUE = '\033[94m'
  OKGREEN = '\033[92m'
  OKCYAN = '\x1b[96m'
  WARNING = '\033[93m'
  FAIL = '\033[91m'
  ENDC = '\033[0m'
  BOLD = '\033[1m'
  UNDERLINE = '\033[4m'


##############
# UI Dialoag #
##############
def exit(text):
  frame = inspect.stack()[1]
  module = inspect.getmodule(frame[0])
  name = module.__name__
  msg = bcolors.FAIL + bcolors.BOLD + name + bcolors.ENDC \
        + bcolors.FAIL + ": " + text + bcolors.ENDC
  sys.exit(msg)
  
def warning(text):
  from setting import quiet
  if not quiet:
    msg = bcolors.WARNING + text + bcolors.ENDC
    print msg

def progress(title, *texts):
  from setting import quiet
  if not quiet:
    msg = bcolors.OKCYAN + bcolors.BOLD + title+":" + bcolors.ENDC
    print msg,
    for info in texts:
      print info,
    sys.stdout.flush()

def done(*texts):
  from setting import quiet
  if not quiet:
    for info in texts:
      print info,
    print " DONE"

def report(title, *texts):
  from setting import quiet
  if not quiet:
    msg = bcolors.OKCYAN + bcolors.BOLD + title+":" + bcolors.ENDC
    print msg,
    for info in texts:
      print info,
    print ""

def prompt(text):
  from setting import no_warning
  if not no_warning:
    frame = inspect.stack()[1]
    module = inspect.getmodule(frame[0])
    name = module.__name__
  
    msg = bcolors.WARNING + name + ": " + bcolors.ENDC
  
    user_input = raw_input(msg + text + \
                 "\nAny key to confirm, enter to cencel...? ")
    if not user_input:
      exit("... ABORT from " + name)
    else:
      report(name, "continue")

def status(title, *texts):
  from setting import quiet
  if not quiet:
    msg = bcolors.OKBLUE + title+":" + bcolors.ENDC
    print msg,
    for info in texts:
      print info,
    print ""
##### END OF UI Diolog #####

###################################
# Simple text formating functions #
###################################
def delete_next(target, pattern, line_number):
  itr = 0
  matched = False
  for line in fileinput.input(target, inplace=1):
    if pattern in line:
      matched = True
    if matched and itr < line_number and itr > 0:
      itr += 1
    else:
      print line,

def delete(target, pattern, line_number):
  itr = 0
  matched = False
  for line in fileinput.input(target, inplace=1):
    if pattern in line:
      matched = True
    if matched and itr < line_number:
      itr += 1
    else:
      print line,

def insert(target, pattern, text):
  for line in fileinput.input(target, inplace=1):
    print line,
    if pattern in line:
      print text

def containing(target, pattern):
  result = False
  with open(target,"r") as ftarget:
    for line in ftarget:
      if pattern in line:
        result = True
  return result
##### END OF text formating #####

#def have_bond(obmol, type_a, type_b):
#  result = False
#  na1 = n2Z(type_a)
#  nb1 = n2Z(type_b)
#  na2 = n2Z(type_b)
#  nb2 = n2Z(type_a)
#  #bond = ob.OBBond()
#  #atom_a = ob.OBAtom()
#  #atom_b = ob.OBAtom()
#  for i in range(obmol.NumBonds()):
#    bond = obmol.GetBond(i)
#    atom_a = bond.GetBeginAtom()
#    atom_b = bond.GetBeginAtom()
#    za = atom_a.GetAtomicNum()
#    zb = atom_b.GetAtomicNum()
#    if (za == na1 and zb == nb1) or (za == na2 and zb == nb2):
#    #print type_a + "-" + type_b + "<=>" + str(za) + "-" + str(zb)
#      result = True
#
#  return result
#  del bond, atom_a, atom_b

#def qt2ob(qtmol):
#  mol = ob.OBMol()
#  for atom in xrange(qtmol.N):
#    new_atom = mol.NewAtom()
#    new_atom.SetAtomicNum(qtmol.Z[atom])
#    new_atom.SetVector(qtmol.R[atom][0], 
#                       qtmol.R[atom][1], 
#                       qtmol.R[atom][2])
#  mol.ConnectTheDots()
#  return mol
#  #del qtmol
#
#def ob2qt(obmol):
#  mol = qg.Molecule()
#  mol.N = obmol.NumAtoms()
#  atom = ob.OBAtom()
#  bond = ob.OBBond()
#
#  _Z = []
#  _coordX = []
#  _coordY = []
#  _coordZ = []
#
#  for i in range(1, obmol.NumAtoms()+1):
#    atom = obmol.GetAtom(i)
#    _Z.append(atom.GetAtomicNum())
#    _coordX.append(atom.GetVector().GetX())
#    _coordY.append(atom.GetVector().GetY())
#    _coordZ.append(atom.GetVector().GetZ())
#  for i in range(obmol.NumBonds()):
#    bond = obmol.GetBond(0)
#    #if bond is not None:
#    obmol.DeleteBond(bond)
#  for i in range(1, obmol.NumAtoms()+1):
#    atom = obmol.GetAtom(1)
#    obmol.DeleteAtom(atom)
#  del obmol
#
#  mol.Z = np.array(_Z)
#  mol.R = np.vstack([_coordX, _coordY, _coordZ]).T
#
#  return mol
#
#  del _coordX, _coordY, coordZ, _Z, mol

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
