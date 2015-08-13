import qctoolkit as qtk
import optimizer as opt
import random
import numpy as np
from time import sleep
import multiprocessing as mp

class MonteCarlo(opt.Optimizer, opt.Temperature):
  def __init__(self, 
               penalty, 
               penalty_input, 
               inpgen,
               **kwargs):
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
  def push(self, penalty, coord):
    self.penalty.append(penalty)
    self.current_penalty = penalty
    qtk.report("MonteCarlo", "step:%d T:%f penalty:%f " % \
               (self.step, self.T, penalty) + \
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
      return np.exp(-abs(dE)/float(self.T))
    if len(args) == 1:
      new_coord = args[0]
    else:
      new_coord = self.getInput()
    out = self.evaluate(new_coord, self.penalty_input)
    if len(self.coord) > 0:
      # accept if out < current_penalty
      # otherwise go through MonteCarlo cycle
      if out >= self.current_penalty:
        rand = random.uniform(0, 1)
        diff = out - self.current_penalty
        boltz = boltzmann(diff)
        qtk.report("MonteCarlo", 
                   "accept worse results?",
                   "rand:%.4f boltz:%.4f dE:%.4f itr:%i"\
                   % (rand, boltz, diff, self.sample_itr))
        # rejection from MonteCarlo cycle, froceed recursion
        if rand > boltz:
          qtk.report("MonteCarlo", "new move rejected", 
                     color='yellow')
          new_coord = self.getInput()
          self.sample_itr += 1 # record MonteCarlo iteration
          if self.parallel == 1 or self.step == 1:
            out, tmp = self.sample()
          # only first finished thread put to queue
          elif self.qout.empty():
            try:
              out, tmp = self.sample()
            # error produced for NoneType return
            except TypeError:
              qtk.report("MonteCarlo", 
                         "job done from another thread")
          # others return None
          else: 
            out = None
            new_coord = None
    if self.parallel == 1 or self.step == 1:
      return out, new_coord
    elif self.qout.empty() and type(out) != None:
      self.qout.put([out, new_coord])

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # parallel process communicate through queue
  def parallel_listener(self,result_queue):
    qtk.progress("MonteCarlo", "waiting results from",
                 self.parallel, "threads")
    # check result every 0.01 second
    while result_queue.empty():
      sleep(0.01)
    result = result_queue.get()
    qtk.done()
    return result[0], result[1]

  def run(self):
    self.step = 1
    init_coord = self.getInput()
    self.current_penalty = self.sample(init_coord)
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
        new_penalty, new_coord = self.parallel_listener(self.qout)
      # serial routine
      else:
        new_penalty, new_coord = self.sample()
      self.push(new_penalty, new_coord)
      self.current_penalty = new_penalty

      
      
