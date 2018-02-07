import qctoolkit as qtk
import optimizer as opt
import random
import numpy as np
import multiprocessing as mp
import operator
import copy

class GeneticOptimizer(opt.Optimizer):
  def __init__(self, penalty, target_input, inpgen, 
               mating_function, pop_size, 
               mutation_rate=0.05, **kwargs):

    opt.Optimizer.__init__(self, penalty, target_input, 
                           inpgen, **kwargs)

    assert pop_size > 2
    self.pop_size = pop_size
    self.mating_function = mating_function
    self.mutation_rate = mutation_rate
    self.old_list = []

  def get_pop(self):
    size = self.pop_size
    if self.mode == 'minimize':
      order = 'ascent'
    elif self.mode == 'maximize':
      order = 'descent'
    else:
      order = 'ascent'
      qtk.warning("mode %s not reconized, set to minimize" % self.mode)
    old_list_db = self.log.list(order=order, has_data=True)[:size]
    old_list = [eval(q.content) for q in old_list_db]
    pop_list = []
    #for i in range(self.threads):
    while len(pop_list) < self.threads:
      if len(old_list) > 2:
        parent1, parent2 = random.sample(old_list, 2)
        pop = self.mating_function(parent1, parent2, self.mutation_rate)
        if not self.repeated(pop):
          pop_list.append(pop)
      else:
        pop = self.getInput()
        if not self.repeated(pop):
          pop_list.append(pop)
    return pop_list

  def fitness(self, pop_list):
    if self.threads == 1:
      fit = []
      for coord in pop_list:
        out = self.evaluate(coord, self.target_input)
        fit.append(out)
      output = np.array(fit).T.tolist()
      return output
    else:
      job = []
      for coord in pop_list:
        job.append([coord, self.target_input])
      out = qtk.parallelize(
        self.evaluate, job, n_output=2, threads=self.threads
      )
      return out

  def run(self, report_step=100):
    step = 0
    qtk.setting.quiet = True
    while not self.converged() and step < self.max_step:
      if report_step:
        if step % report_step == 0:
          qtk.setting.quiet = False
          qtk.progress('GA', '%d steps' % step)
      pop = self.get_pop()
      qtk.progress("Optimizer", "GE iteration with %d new points" % len(pop))
      self.register(pop)
      fit, info = self.fitness(pop)
      step += 1
      if type(fit) is list:
        self.update(pop, fit, info)
      else:
        self.update(pop, [fit], [info])
      qtk.setting.quiet = True
