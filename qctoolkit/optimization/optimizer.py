import qctoolkit as qtk
import numpy as np

class Temperature(object):
  """
  temperature object for temerature dependent optimizer
  e.g. MonteCarlo, basinhopping, simulated annealing...
  """
  def __init__(self, T, **kwargs):
    self.T = T
    if 'annealing' in kwargs:
      self.annealing = kwargs['annealing']
    else:
      self.annealing = 'exp'
    if self.annealing and 'k_factor' not in kwargs:
      if self.annealing == 'exp':
        self.k = 1
      elif self.annealing == 'fast':
        self.k = 1.1
      elif self.annealing == 'boltz':
        self.k = 2.75
      else: qtk.exit(self.annealing, " is not implemented")
    elif self.annealing:
      self.k = kwargs['k_factor']

  # gradually decrease T to zero, simulated annealing
  def decrease_T(self):
    if self.annealing:
      if self.annealing == 'exp':
        self.T = self.T * (0.95)**self.k
      elif self.annealing == 'fast':
        self.T = self.T/self.k
      else:
        self.T = self.T/np.log(self.k)

  # gradually increase T, for minima/basin hopping
  def increase_T(self):
    pass

class Optimizer(object):
  """
  general object for optimization
  penalty_funciton: function to be minimized to 0
  input_generator: function to generate new input
  """
  def __init__(self, 
               penalty_function, 
               penalty_input,
               input_generator,
               **kwargs):

    ################################
    # setting default setup values #
    ################################
    if 'cutoff' in kwargs:
      self.cutoff = kwargs['cutoff']
    else:
      self.cutoff = 1E-5

    if 'max_step' in kwargs:
      self.max_step = kwargs['max_step']
    else:
      self.max_step = 1000

    # set abs(penalty)**power
    # if power==0, no modification is applied
    if 'power' in kwargs:
      self.power = kwargs['power']
    else:
      self.power = 0

    if 'log_file' in kwargs:
      self.log = kwargs['log_file']
    else:
      self.log = 'optimization.log'

    if 'average_length' in kwargs:
      self.avg_len = kwargs['average_length']
    else:
      self.avg_len = 100

    if 'dump' in kwargs:
      self.dump_len = kwargs['dump']
    else:
      self.dump_len = 30
    ##### end of default values setup #####

    # open logfile, close when converged
    self.logfile = open(self.log, 'w', 0)
    self.logfile.flush()
    # result lists
    self.penalty = []
    self.coord = []
    self.step = 0
    # parsing variables
    self.penalty_function = penalty_function
    self.penalty_input = penalty_input
    self.input_generator = input_generator

  ################################################
  # GENERAL UTILITY FUNCTIONS OF EVERY OPTIMIZER #
  ################################################

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # each function implement run differently
  def run(self):
    raise NotImplementedError("Please Implement run method")

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # flush penalty/coord list to log file
  def dump(self):
    size = self.dump_len
    output = list(zip(self.penalty[:size], self.coord[:size]))
    self.penalty = self.penalty[size:]
    self.coord = self.coord[size:]
    for line in output:
      print >> self.logfile, "% 10.6E %s" % (line[0], line[1])
      self.logfile.flush()

  # !!!!!!!!!!!!!!!!!!!!!!!!!!
  # update penalty/coord lists
  def push(self, penalty, coord):
    self.step += 1
    self.penalty.append(penalty)
    self.coord.append(coord)
    if self.step > 0 \
     and self.step % self.dump_len ==0\
     and len(self.coord) > 3*self.dump_len:
      self.dump()

  # !!!!!!!!!!!!!!!!!!!!!!!!!!
  # convergence/max_step check
  def converged(self):
    size = self.avg_len
    length = min(len(self.penalty), size)
    if length > 1:
      if sum(self.penalty[-length:])/float(length)\
           < self.cutoff\
         or self.step >= self.max_step:
        self.logfile.close()
        return True
    else: return False

  # !!!!!!!!!!!!!!!!!!!!!!
  # evaluate penalty value
  def evaluate(self, new_coord, penalty_input):
    if type(penalty_input[-1]) == dict:
      evl_args = penalty_input[:-1]
      evl_kwgs = penalty_input[-1]
      out = self.penalty_function(\
              new_coord, *evl_args, **evl_kwgs)
    else:
      out = self.penalty_function(new_coord, *evl_args)
    if self.power == 0:
      return out
    else:
      return abs(out)**self.power

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # generate input for penalty function by input_generator
  def getInput(self, *inp_args, **inp_kwargs):
    return self.input_generator(*inp_args, **inp_kwargs)
