import qctoolkit as qtk
import optimizer as opt
import random
import numpy as np
import multiprocessing as mp
import operator

class GeneticOptimizer(opt.Optimizer):
  def __init__(self, penalty, penalty_input, inpgen, 
               mating_function, pop_size, **kwargs):
    # redefine some functions are necessary!
    opt.Optimizer.__init__(self, penalty, penalty_input, 
                           inpgen, **kwargs)
    self.pop_size = pop_size
    self.mating_function = mating_function
    self.mutation_rate = 0.05

    if not pop_size/2 > 2:
      print pop_size
      qtk.exit("population too small")

  def get_pop(self, size):
    pop_list = []
    for i in range(size):
      pop_list.append(self.getInput())
    return pop_list

  def fitness(self, pop_list):
    fit = []
    # serial implementation
    if self.parallel == 1:
      for coord in pop_list:
        fit.append(self.evaluate(coord, self.penalty_input))
    else:

      def run_job(qinp, qout):
        while not qinp.empty():
          inpd = qinp.get()
          inp = inpd[0]
          ind = inpd[1]
          out = self.evaluate(inp, self.penalty_input)
          qout.put([out, ind])

      qout = mp.Queue()
      qinp = mp.Queue()
      for i in range(len(pop_list)):
        inp = [pop_list[i], i]
        qinp.put(inp)
      proc = []
      for i in range(self.parallel):
        p = mp.Process(target = run_job, args=(qinp,qout))
        p.start()
        proc.append(p)
      for p in proc:
        p.join()

      for i in range(len(pop_list)):
        fit.append(qout.get())
      fit = list(np.array(sorted(
              fit, key=operator.itemgetter(1)))[:,0])
    return fit

  def pop_sort(self):
    pop_sorted = [[fit, pop] for fit, pop in\
                sorted(zip(self.fit, self.pop))]
    self.pop = [pop[1] for pop in list(pop_sorted)]
    self.fit = [pop[0] for pop in list(pop_sorted)]

  def merge(self, new_pop, new_fit):
    self.pop.extend(new_pop)
    self.fit.extend(new_fit)

  def get_newpop(self):
    self.pop = self.pop[:self.pop_size/2]
    self.fit = self.fit[:self.pop_size/2]
    index = range(len(self.pop))
    new_pop = []
    for i in range(self.pop_size - len(index)):
      ind = random.sample(index, 2)
      parent1 = self.pop[ind[0]]
      parent2 = self.pop[ind[1]]
      new_pop.append(self.mating_function(\
                     parent1, parent2, self.mutation_rate\
                    ))
    return new_pop
    
  def run(self):
    self.pop = self.get_pop(self.pop_size)
    self.fit = self.fitness(self.pop)
    self.pop_sort()
    self.push(self.fit[0], self.pop[0])
    while not self.converged():
      new_pop = self.get_newpop()
      new_fit = self.fitness(new_pop)
      self.merge(new_pop, new_fit)
      self.pop_sort()
      self.push(self.fit[0], self.pop[0])
      
    


  # construct initial population

  # evaluate fitness (penalty) function 

  # while converge/step loop 
  #  generate new population by crossing
  #
  #  evaluate new fitness function
  #
  #  select best from new/old population -> old population
  #
  # 

    
