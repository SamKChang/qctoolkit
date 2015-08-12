import qctoolkit as qtk
import optimizer as opt
import random
import numpy as np

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
    else:
      _T = 1
    opt.Temperature.__init__(self,_T, **kwargs)

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # override push for simulated annealing
  def push(self, penalty, coord):
    self.step += 1
    self.penalty.append(penalty)
    self.current_penalty = penalty
    if self.annealing:
      self.decrease_T()
      self.coord.append([coord, {'T':self.T}])
    else:
      self.coord.append(coord)
    if self.step > 0 and self.step % self.dump_len ==0:
      qtk.report("MonteCarlo", "step:%d T:%f penalty:%f\n" % \
                 (self.step, self.T, penalty) + \
                 "coord:%s" % coord)
      self.dump()

  def boltzmann(self, dE):
    return np.exp(-abs(dE)/float(self.T))

  def sample(self, new_coord):
    out = self.evaluate(new_coord, self.penalty_input)
    if len(self.coord) > 0:
      if out < self.current_penalty:
        self.push(out, new_coord)
        return out
      else:  
        rand = random.uniform(0, 1)
        diff = out - self.current_penalty
        if rand <= self.boltzmann(diff):
          self.push(out, new_coord)
          return out
        else:
          new_coord = self.getInput()
          out = self.sample(new_coord)
    else:
      self.push(out, new_coord)
      return out
    
  def run(self):
    init_coord = self.getInput()
    self.current_penalty = self.sample(init_coord)
    while not self.converged():
      coord = self.getInput()
      self.sample(coord)
      
      
