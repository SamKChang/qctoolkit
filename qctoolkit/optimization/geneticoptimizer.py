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
    old_list_db = self.log.list(order='ascent')[:size]
    old_list = [eval(q.content) for q in old_list_db]
    pop_list = []
    for i in range(self.threads):
      if len(old_list) > 2:
        parent1, parent2 = random.sample(old_list, 2)
        pop_list.append(self.mating_function(\
          parent1, parent2, self.mutation_rate\
        ))
      else:
        pop_list.append(self.getInput())
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

  def run(self):
    while not self.converged():
      pop = self.get_pop()
      fit, info = self.fitness(pop)
      if type(fit) is list:
        for i in range(len(fit)):
          content = str(pop[i])
          self.log.push(content, fit[i], info[i])
      else:
        content = str(pop)
        self.log.push(content, fit, info)
