import sys, os, glob, re, copy
import qmio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing as mp
import collections as cl

class QMResults(object):
  def __init__(self, pathPattern):
    self.name = ''
    self.unit = 'hartree'
    self.path = re.sub('/$', '', pathPattern[0])
    self.pattern = pathPattern[1]
    self.program = pathPattern[2]
    if len(pathPattern)>3:
      flag = re.compile('(True|False)')
      if re.match(flag, pathPattern[3]):
        self.parallel = pathPattern[3]
      else:
        sys.exit("Error from analysis.py->QMResults: '" +\
                 pathPattern[3] + "' is not a boolean variable "+\
                 "for parallel setup. Please choose either " +\
                 "'True' or 'False'. 'False is default'"
                )
    else:
      self.parallel = False
    self.out_dir = glob.glob(self.path + "/" + self.pattern)
    #self.out_dir = next(os.walk(path))[1]
    #self.out_dir.sort()
    self.results = {}
    self.Et = {}
    self.data = np.atleast_2d(np.array([]))
    if self.path+'inp' in self.out_dir: 
      self.out_dir.remove(self.path + 'inp')

    if self.parallel:
      # parallel with minimal fork
      # queue is necessary for shared data
      def read_out_dir(out_path, Et_queue):
        for folder in out_path:
          for out in glob.glob(folder + '/*.out'):
            # save data of each process in the queue
            data = qmio.QMOut(out, self.program)
            Et_queue.put(
              {re.sub(self.path + "/", "", out) : 
               np.array([data.Et, data.SCFStep])}
            )
  
      def chunks(l, n):
        for i in xrange(0, len(l), n):
          yield l[i:i+n]
  
      job_chunk = list(chunks(
                    self.out_dir, 
                    len(self.out_dir)/mp.cpu_count()
                  ))
  
      jobs = []
      queue = mp.Queue()
      for result in job_chunk:
        p = mp.Process(target=read_out_dir,\
          args=(result,queue))
        jobs.append(p)
        # start multiprocess
        p.start()
      for process in job_chunk:
        for result in process:
          # append data from queue to hash
          self.results.update(queue.get())
      for process in jobs:
        # wait for each process to finish
        process.join()
      
    else:
      # Serial verision
      def read_out_dir(out_path):
        for out in glob.glob(out_path + '/*.out'):
          data = qmio.QMOut(out, self.program)
          self.results.update({
            re.sub(self.path + "/", "", out) :
            np.array([data.Et, data.SCFStep])
          })
          #data.Et
  
      # Serial verision
      for results in self.out_dir:
        read_out_dir(results)

    # reture sorted hash as default
    self.results = cl.OrderedDict(sorted(self.results.iteritems()))
    self.Et = cl.OrderedDict(sorted({
                k : v[0] for k, v in self.results.iteritems()
              }.iteritems()))

  def sum(self):
    return sum(self.Et.itervalues())

  def ev(self):
    if re.match('hartree', self.unit):
      self.Et = { k : float(v) * 27.21138505\
                  for k,v in self.Et.iteritems()}
      self.unit = 'ev'
    elif re.match('kcal', self.unit):
      self.Et = { k : float(v) * 0.0433634\
                  for k,v in self.Et.iteritems()}
      self.unit = 'ev'

  def kcal(self):
    if re.match('hartree', self.unit):
      self.Et = { k : float(v) * 627.509469\
                  for k,v in self.Et.iteritems()}
      self.unit = 'kcal'
    elif re.match('ev', self.unit):
      self.Et = { k : float(v) * 23.061\
                  for k,v in self.Et.iteritems()}
      self.unit = 'kcal'

  def setStep(self, step):
    for key in self.results.keys():
      if self.results[key][1] > step:
        del self.Et[key]
        del self.results[key]

  def rmKey(self, pattern):
    p1 = re.compile(pattern)
    self.results = {
      re.sub(p1, '', k) : v for k, v in self.results.iteritems()
    }
    self.Et = {
      re.sub(p1, '', k) : v for k, v in self.Et.iteritems()
    }
    
  def exKey(self, pattern):
    p1 = re.compile(pattern)
    self.results = {
      p1.match(k).group(1) : v for k, v in self.results.iteritems()
    }
    self.Et = {
      p1.match(k).group(1) : v for k, v in self.Et.iteritems()
    }

  def subtract(self, other_result):
    out = copy.deepcopy(self)
    for key in self.Et:
      if key in other_result.Et:
        out.Et[key] = self.Et[key] - other_result.Et[key]
        out.results[key][0] = \
          self.results[key][0] \
         -other_result.results[key][0]
        out.results[key][1] = max(
                                out.results[key][1],
                                other_result.results[key][1])
      else:
        del out.Et[key]
    return out

  def subtract_constant(self, const):
    out = copy.deepcopy(self)
    out.Et = { k : float(v)-float(const)\
                for k, v in self.Et.iteritems()}
    for key in out.results:
      out.results[key] = np.array([
        self.results[key][0] - float(const),
        self.results[key][1]
      ])
    return out

  def npdata(self):
    try:
      out = np.transpose(np.vstack([
        np.fromiter(self.Et.iterkeys(), dtype=float),
        np.fromiter(self.Et.itervalues(), dtype=float)
      ]))
      self.data = out[out[:,0].argsort()]
    except ValueError:
      #sys.exit("value error, not numeric")
      i=1
      results = []
      for v in self.Et.itervalues():
        results.append([i, v])
        i += 1
      self.data = np.array(results)
    return self.data

  # print total energy in plan text form
  def write(self, value, out_file):
    out = sys.stdout if re.match('stdout', out_file) else\
          open(out_file, "w")
    if re.match(re.compile('^E$'), value):
      for key in self.Et:
        print >> out, "%s\t%s" % (key, self.Et[key])
    else:
      print "other"
      for key in self.results:
        print >> out, "%s\t%s" % (key, self.results[key])

class ScatterPlot(QMResults):
  def __init__(self, QMOut_pred, QMOut_true):
    self.pred = QMOut_pred
    self.true = QMOut_true
    self.data = []
    self.list = []
    match = False
    for key in self.pred.Et:
      if key in self.true.Et:
        self.list.append(key)
        new_point = [self.pred.Et[key], self.true.Et[key]]
        self.data.append(new_point)
        match = True
    self.data = np.array([self.data])[0]
    if match:
      x = self.data[:,0]
      y = self.data[:,1]
      self.fit = np.polyfit(x, y, 1)
      self.y_fit = self.fit[0]*x + self.fit[1]
      diff = abs(y-self.y_fit)
      self.MAE = sum(diff)/len(diff)
      self.RMSE = np.sqrt(sum(np.square(diff))/len(diff))
      #print self.MAE, self.RMSE, self.list
    else:
      sys.exit("ERROR from analysis.py->ScatterPlot: "+\
               "data keys not matched"
              )

  def plot(self):
    plot_data = np.hstack([
      self.data, np.transpose(np.atleast_2d(self.y_fit))
    ])
    plot_data = plot_data[plot_data[:,0].argsort()]
    x = plot_data[:,0]
    y1 = plot_data[:,1]
    y2 = plot_data[:,2]

    self.fig = plt.figure(figsize=(9, 8))
    ax = self.fig.add_subplot(1,1,1, adjustable='box', aspect=1)
    ax.plot(x, y1, 'ko', 
      markersize = 10,
      markeredgewidth = 1,
      markerfacecolor = 'none'
    )
    ax.plot(x, y2, 'r--',
      linewidth = 1.5
    )
    ax.set_xlabel(self.pred.name, fontsize=20)
    ax.set_ylabel(self.true.name, fontsize=20)
    ax.tick_params(labelsize=15)

  def show(self):
    plt.show()

  def save(self, out_file):
    self.fig.savefig(out_file)
    

class DensityPlot(QMResults):
  def __init__(self, path):
    self.out_dir = QMResults(path)
