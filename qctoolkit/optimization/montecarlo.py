import qctoolkit as qtk
import optimizer as opt
import random
import numpy as np
from time import sleep
import multiprocessing as mp

class MonteCarlo(opt.Optimizer, opt.Temperature):
  """
  penalty: function to optimize to 0
  penalty_input: list of all input arguments of penalty function
                 exept coordinate, which will be randomly sampled
                 and put as the 'FIRST' argument of the 
                 penalty function
  inpgen: input generator, for MonteCarlo, it takes no input
  """
  def __init__(self, penalty, penalty_input, inpgen, **kwargs):
    opt.Optimizer.__init__(self, penalty, penalty_input, 
                           inpgen, **kwargs)
    if 'T' in kwargs:
      _T = kwargs['T']
      del kwargs['T']
    else:
      _T = 1
    opt.Temperature.__init__(self,_T, **kwargs)
    self.sample_itr = 0

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # override push for simulated annealing
  def push(self, penalty, result, coord):
    self.penalty.append(penalty)
    self.result.append(result)
    self.current_penalty = penalty
    qtk.report("MonteCarlo", "step:%d T:%f penalty:%f result:%f "%\
               (self.step, self.T, penalty, result) + \
               "coord:%s" % coord, color='green')
    if self.annealing:
      self.decrease_T()
      self.coord.append([coord, {'T':"%.3E" % self.T}])
    else:
      self.coord.append(coord)
    self.write_log()
    if self.step > 0 \
    and self.step % self.dump_len ==0\
    and len(self.coord) > 3*self.dump_len:
      self.dump()


  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # MonteCarlo sampling function, MAIN ROUTINE
  def sample(self, *args):
    def boltzmann(dE):
      if self.T > 0:
        return np.exp(-abs(dE)/float(self.T))
      else:
        return 0
    if len(args) == 1:
      new_coord = args[0]
    else:
      new_coord = self.getInput()
    penalty, out = self.evaluate(new_coord, self.penalty_input)
    if len(self.coord) > 0:
      # accept if penalty < current_penalty
      # otherwise go through MonteCarlo cycle
      if penalty >= self.current_penalty:
        rand = random.uniform(0, 1)
        diff = penalty - self.current_penalty
        boltz = boltzmann(diff)
        qtk.report("MonteCarlo", 
                   "accept worse results?",
                   "rand:%.4f boltz:%.4f dE:%.4f itr:%i"\
                   % (rand, boltz, diff, self.sample_itr))
        # rejection from MonteCarlo cycle, froceed recursion
        if rand > boltz:
          qtk.report("MonteCarlo", "new move rejected", 
                     color='yellow')
          new_coord = self.getInput() # generate new random inp
          self.sample_itr += 1 # record MonteCarlo iteration
          if self.parallel == 1 or self.step == 1:
            # iterative run for serial job
            penalty, out, _ = self.sample()
          # only first finished thread put to queue
          elif self.qout.empty():
            # iterative run for parallel case
            try:
              penalty, out, _ = self.sample()
              print penalty, out
            # error produced for NoneType return
            except TypeError:
              qtk.report("MonteCarlo", 
                         "job done from another thread")
          # others return None
          else: 
            out = None
            penalty = None
            new_coord = None
    # serial return
    if self.parallel == 1 or self.step == 1:
      return penalty, out, new_coord
    # parallel case, put to queue instead of return
    elif self.qout.empty() and type(penalty) != None:
      self.qout.put([penalty, out, new_coord])

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # parallel process communicate through queue
  # regularly check ouput queue for result
  # return results from first accepted thread to master 
  # 
  # Communication mechanism: Using self.qout
  #  accepted result is put to queue
  #  only first accepted thread is put to queue
  #  queue is always empty, waiting for accepted result
  #  parallel_listener moniters queue and return result
  def parallel_listener(self,result_queue):
    qtk.progress("MonteCarlo", "waiting results from",
                 self.parallel, "threads")
    # check result every 0.01 second, HARD CODED
    while result_queue.empty():
      sleep(0.01)
    result = result_queue.get() # collect results
    qtk.done() # indicate job done
    # return results to master
    return result[0], result[1], result[2] 

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # implementation of run method from Optimizer object
  def run(self):
    # initial run
    self.step = 1
    init_coord = self.getInput()
    self.current_penalty = self.sample(init_coord)
    # loop over convergence/max_step
    while not self.converged():
      self.step += 1
      self.sample_itr = 0
      # parallel routine
      if self.parallel > 1:
        self.qout = mp.Queue()
        procs = []
        for _ in range(self.parallel):
          p = mp.Process(target=self.sample)
          p.start()
          procs.append(p)
        for p in procs:
          p.join()
        # collect result from single finished thread
        new_penalty, new_out, new_coord = \
          self.parallel_listener(self.qout)
      # serial routine
      else:
        new_penalty, new_out, new_coord = self.sample()
      # update result
      self.push(new_penalty, new_out, new_coord)
      self.current_penalty = new_penalty
