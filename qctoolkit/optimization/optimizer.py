import qctoolkit as qtk
import numpy as np
import os

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
  target_function: function to be minimized to 0
  input_generator: function to generate new input
  """
  def __init__(self, 
               target_function, 
               target_input,
               input_generator,
               cutoff=1E-5,
               power=1,
               threads=1,
               target=0,
               converge_length=20,
               distributed=False,
               log='optimization.db',
               new_run=True,
               **kwargs):

    ################################
    # setting default setup values #
    ################################
    self.cutoff = cutoff
    self.power = power
    if os.path.exists(log):
      if new_run:
        os.remove(log)
    self.log = qtk.Logger(log)
    self.threads = threads
    self.target = target
    self.converge_length = converge_length
    self.conv_itr = 0
    self.distributed=distributed
    ##### end of default values setup #####

    # result lists
    self.result = []
    self.penalty = []
    self.coord = []
    self.step = 0
    # parsing variables
    self.target_function = target_function
    self.target_input = target_input
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
    content = str(self.coord[-1]) + ', ' + str(self.result[-1])
    self.log.push(content, self.penalty[-1])

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

  # !!!!!!!!!!!!!!!!!!!!!!!!!!
  # convergence/max_step check
  def converged(self):
    # first step is set to NOT converge
    penalty = [q.data for q in self.log.list()]
    if len(penalty)>1:
      # if get to target
      if penalty[-1] < self.cutoff:
        #self.logfile.close()
        return True
      # if stuck at best solution
      elif penalty[-1] == penalty[-2]:
        if self.conv_itr > self.converge_length:
          #self.logfile.close()
          return True
        else:
          self.conv_itr = self.conv_itr + 1
          return False
      # if reach max step
      elif self.step >= self.max_step:
        qtk.report("Optimizer", "max_step reached stopping",
                   color='red')
        #self.logfile.close()
        return True
      else: 
        self.conv_itr = 0
        return False
    else: return False

  # !!!!!!!!!!!!!!!!!!!!!!
  # evaluate penalty value
  def evaluate(self, new_coord, target_input):
    if type(target_input[-1]) == dict:
      evl_args = target_input[:-1]
      evl_kwgs = target_input[-1]
      out = self.target_function(\
              new_coord, *evl_args, **evl_kwgs)
    else:
      out = self.target_function(new_coord, *evl_args)
    if type(out) is tuple or type(out) is list:
      out_info = out[1:]
      out = out[0]
    else:
      out_info = out
    return abs(out - self.target)**self.power, str(out_info)

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # generate input for penalty function by input_generator
  def getInput(self, *inp_args, **inp_kwargs):
    return self.input_generator(*inp_args, **inp_kwargs)
