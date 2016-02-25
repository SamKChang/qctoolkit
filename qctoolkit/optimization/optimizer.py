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
      self.power = 1

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

    if 'parallel' in kwargs:
      self.parallel = kwargs['parallel']
    else:
      self.parallel = 1

    if 'target' in kwargs:
      self.target = kwargs['target']
    else:
      self.target = 0

    if 'converge_length' in kwargs:
      self.converge_length = kwargs['converge_length']
    else:
      self.converge_length = 20
    
    self.conv_itr = 0
    self.log_step = 1


    ##### end of default values setup #####

    # open logfile, close when converged
    self.logfile = open(self.log, 'w', 0)
    self.logfile.flush()
    # result lists
    self.result = []
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
  def write_log(self):
    if self.step % self.log_step ==0:
      line = [self.result[-1], self.penalty[-1], self.coord[-1]]
      print >> self.logfile, "% 10.6E %s" % (line[0], line[2])
      self.logfile.flush()

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
  def push(self, penalty, result, coord):
    self.step += 1
    self.penalty.append(penalty)
    self.result.append(result)
    self.coord.append(coord)
    qtk.report("Optimizer", 
               "result:%.4E penalty:%.4E coord%s itr:%d"\
               % (result, penalty, coord, self.step))
    self.write_log()
    if self.step > 0 \
     and self.step % self.dump_len ==0\
     and len(self.coord) > 3*self.dump_len:
      self.dump()

  # !!!!!!!!!!!!!!!!!!!!!!!!!!
  # convergence/max_step check
  def converged(self):
    # first step is set to NOT converge
    if len(self.penalty)>1:
      # if get to target
      if self.penalty[-1] < self.cutoff:
        self.logfile.close()
        return True
      # if stuck at best solution
      elif self.penalty[-1] == self.penalty[-2]:
        if self.conv_itr > self.converge_length:
          self.logfile.close()
          return True
        else:
          self.conv_itr = self.conv_itr + 1
          return False
      # if reach max step
      elif self.step >= self.max_step:
        qtk.report("Optimizer", "max_step reached stopping",
                   color='red')
        self.logfile.close()
        return True
      else: 
        self.conv_itr = 0
        return False
    else: return False

    # average over length...
#    size = self.avg_len
#    length = min(len(self.penalty), size)
#    if length > 1:
#      if sum(self.penalty[-length:])/float(length) < self.cutoff\
#      or self.step >= self.max_step:
#        if self.step >= self.max_step:
#          qtk.report("Optimizer", "not max_step reached stop",
#                     color='red')
#        self.logfile.close()
#        return True
#    else: return False

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
    return abs(out - self.target)**self.power, out

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # generate input for penalty function by input_generator
  def getInput(self, *inp_args, **inp_kwargs):
    return self.input_generator(*inp_args, **inp_kwargs)
