from math import pi ,sin, cos
import molecule as qg
#import openbabel as ob
import numpy as np
#import gc
import fileinput
import sys, re, copy
import inspect
import multiprocessing as mp
import operator
from compiler.ast import flatten
import qctoolkit.data.elements as qel

def convE(source, units):

  def returnError(ioStr, unitStr):
    msg = 'supported units are:\n'
    for key in Eh.iterkeys():
      msg = msg + key + '\n'
    report(msg, color=None)
    exit(ioStr + " unit: " + unitStr + " is not reconized")

  EhKey = {
    'eh': 'Eh',
    'ev': 'eV',
    'kcal/mol': 'kcal/mol',
    'cminv': 'cmInv',
    'k': 'K',
    'j': 'J',
    'kj/mol': 'kJ/mol',
    'joule': 'J',
    'Joule': 'J',
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

  unit = units.split('-')
  if unit[0].lower() != 'hartree' and unit[0].lower() != 'eh':
    if unit[0].lower() in EhKey:
      unit0 = EhKey[unit[0].lower()]
      source = source / Eh[unit0]
    else: returnError('input', unit[0])
  if unit[1].lower() not in EhKey: 
    returnError('output', unit[1])
  else:
    unit1 = EhKey[unit[1].lower()]
  return source * Eh[unit1]


def fileStrip(path):
  new_path = re.sub('.*/', '', path)
  return new_path

def imported(module):
  try:
    __import__(module)
  except ImportError:
    return False
  else:
    return True

def Structure(input_data, **kwargs):
  if type(input_data) is not qg.Molecule:
    try:
      return qg.Molecule(input_data)
    except:
      pass
  else:
    return input_data

def pathStrip(path):
  new_path = re.sub('//*','/',path)
  new_path = re.sub(r'([^\.])\.\/',r"\1",new_path)
  return new_path

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

def parallelize(target_function, 
                input_list, 
                **kwargs):
  import setting
 
  if 'threads' in kwargs:
    threads = kwargs['threads']
  else:
    threads = setting.cpu_count
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
    for i in range(0, len(_list), _size):
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
  raise RuntimeError(msg)
  sys.exit(msg)
  
def warning(text):
  from setting import quiet
  if not quiet:
    msg = bcolors.WARNING + text + bcolors.ENDC
    print msg
  sys.stdout.flush()

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
  sys.stdout.flush()

def report(title, *texts, **kwargs):
  from setting import quiet
  if not quiet:
    if 'color' in kwargs:
      color = kwargs['color']
    else:
      color = 'cyan'
    tle = bcolors.ENDC
    if color == 'cyan':
      msghead = bcolors.OKCYAN + bcolors.BOLD
    elif color == 'blue':
      msghead = bcolors.OKBLUE + bcolors.BOLD
    elif color == 'green':
      msghead = bcolors.OKGREEN + bcolors.BOLD
    elif color == 'yellow':
      msghead = bcolors.WARNING + bcolors.BOLD
    elif color == 'red':
      msghead = bcolors.FAIL + bcolors.BOLD
    else:
      msghead = ''
      tle = ''

    msg = msghead + title + ":" + tle
    print msg,
    for info in texts:
      print info,
    print ""
  sys.stdout.flush()

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
  sys.stdout.flush()

def status(title, *texts):
  from setting import quiet
  if not quiet:
    msg = bcolors.OKBLUE + title+":" + bcolors.ENDC
    print msg,
    for info in texts:
      print info,
    print ""
  sys.stdout.flush()
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

def replace(target, pattern, text):
  for line in fileinput.input(target, inplace=1):
    print re.sub(re.compile(pattern), text, line),

def containing(target, pattern):
  result = False
  with open(target,"r") as ftarget:
    for line in ftarget:
      if pattern in line:
        result = True
  return result

def matching(target, pattern):
  result = False
  with open(target,'r') as ftarget:
    for line in ftarget:
      if re.match(re.compile(pattern),line):
        result = True
  return result
##### END OF text formating #####

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
    exit("n2ve: element type " + Zn + " is not defined")

def Z2n(Z):
  if type_list.has_key(Z):
    return type_list[Z]
  else:
    exit("Z2n: atomic number " + str(Z) + " is not defined")
    #return Z
  
def n2Z(Zn):
  if z_list.has_key(Zn):
    return z_list[Zn]
  else:
    exit("n2Z: element type " + str(Zn) + " is not defined")
  
def n2m(Zn):
  if mass_list.has_key(Zn):
    return mass_list[Zn]
  else:
    exit("n2Z: element type " + str(Zn) + " is not defined")

def qAtomName(query):
  if type(query) == str:
    if z_list.has_key(query):
      return str(query)
  elif type(query) == int or type(query) == float:
    if type_list.has_key(int(query)):
      return str(Z2n(query))
  else:
    exit("qAtom: element " + str(Zn) + " is not defined")

def qAtomicNumber(query):
  if type(query) == str:
    if z_list.has_key(query):
      return n2Z(query)
  elif type(query) == int or type(query) == float:
    if type_list.has_key(int(query)):
      return query
  else:
    exit("qAtom: element " + str(Zn) + " is not defined")

def isAtom(query):
  if z_list.has_key(query):
    return True
