import sys, os, glob, re, copy
import qmio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import multiprocessing as mp
import collections as cl
import pandas as pd
import utilities as ut

# Construct dictionary of {file:results} for DataFrame
class QMResults(object):
  def __init__(self, path, pattern, program, **kwargs):
    self.path = re.sub('/$', '', path)
    self.pattern = pattern
    self.program = program
    if 'threads' in kwargs:
      self.threads = int(kwargs['threads'])
    else:
      self.threads = 1

    self.out_dir = glob.glob(self.path + "/" + self.pattern)
    print "Reading output file in " +\
          self.path + "/" + self.pattern
    self.results = {}
    self.data = np.atleast_2d(np.array([]))
    if self.path+'inp' in self.out_dir: 
      self.out_dir.remove(self.path + 'inp')

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
                  len(self.out_dir)/self.threads
                ))

    jobs = []
    itr = 1
    queue = mp.Queue()
    for result in job_chunk:
      p = mp.Process(target=read_out_dir,\
        args=(result,queue))
      jobs.append(p)
      # start multiprocess
      p.start()
      ut.report("reading output", "thead", itr, "started...")
      itr += 1
    itr = 0
    for process in job_chunk:
      ut.report("collecting results", "from thread:", itr)
      for result in process:
        # append data from queue to hash
        self.results.update(queue.get())
      itr += 1
    for process in jobs:
      # wait for each process to finish
      ut.report("waiting", "from thread:", itr)
      process.join()

# pandas DataFrame wrapper
class QMData(pd.DataFrame):
  def __init__(self, arg_path, arg_pattern, arg_prog, **kwargs):
    if 'threads' in kwargs:
      arg_threads = kwargs['threads']
    else:
      arg_threads = 1
    _qr = QMResults(arg_path, 
                    arg_pattern, 
                    arg_prog,
                    threads=arg_threads).results
    _qr = pd.DataFrame(_qr).T
    _qr.index.name = 'file'
    _qr.columns = ['E', 'step']
    super(QMData, self).__init__(_qr)
    #self.E_unit = 'Ha'

  ##########################
  # extract regex of index #
  ##########################
  def extract(self, pattern):
    self.index = self.index.to_series().astype(str)\
                 .str.extract(pattern)

  ####################################
  # mathematical operations to index #
  ####################################
  def index_add(self, factor):
    self.index = self.index.to_series().astype(float) + factor

  def index_sub(self, factor):
    self.index = self.index.to_series().astype(float) - factor

  def index_mul(self, factor):
    self.index = self.index.to_series().astype(float) * factor

  def index_div(self, factor):
    self.index = self.index.to_series().astype(float) / factor


  #######################
  #  OPERATOR OVERLOAD  #
  #######################
  def __add__(self, other):
    _out = copy.deepcopy(self)
    if isinstance(other, QMData):
      _term1 = copy.deepcopy(pd.DataFrame(self))
      _term2 = copy.deepcopy(pd.DataFrame(other))
      _term2.columns = ['E2', 'step2']
      _tmp = pd.concat([_term1, _term2], axis=1)\
               [['step','step2']].max(axis=1)
      _tmp.columns = ['step']
      _out['step'] = _tmp
      _term2.columns = ['E', 'step']
      _term3 = _term1 + _term2
      _out['E'] = _term3['E']
    elif isinstance(other, float) or isinstance(other, int):
      _out['E'] = self['E'] + other
    return _out

  def __sub__(self, other):
    _out = copy.deepcopy(self)
    if isinstance(other, QMData):
      _term1 = copy.deepcopy(pd.DataFrame(self))
      _term2 = copy.deepcopy(pd.DataFrame(other))
      _term2.columns = ['E2', 'step2']
      _tmp = pd.concat([_term1, _term2], axis=1)\
               [['step','step2']].max(axis=1)
      _tmp.columns = ['step']
      _out['step'] = _tmp
      _term2.columns = ['E', 'step']
      _term3 = _term1 - _term2
      _out['E'] = _term3['E']
    elif isinstance(other, float) or isinstance(other, int):
      _out['E'] = self['E'] - other
    return _out

  def __mul__(self, other):
    _out = copy.deepcopy(self)
    if isinstance(other, QMData):
      _term1 = copy.deepcopy(pd.DataFrame(self))
      _term2 = copy.deepcopy(pd.DataFrame(other))
      _term2.columns = ['E2', 'step2']
      _tmp = pd.concat([_term1, _term2], axis=1)\
               [['step','step2']].max(axis=1)
      _tmp.columns = ['step']
      _out['step'] = _tmp
      _term2.columns = ['E', 'step']
      _term3 = _term1 * _term2
      _out['E'] = _term3['E']
    elif isinstance(other, float) or isinstance(other, int):
      _out['E'] = self['E'] * other
    return _out

  def __div__(self, other):
    _out = copy.deepcopy(self)
    if isinstance(other, QMData):
      _term1 = copy.deepcopy(pd.DataFrame(self))
      _term2 = copy.deepcopy(pd.DataFrame(other))
      _term2.columns = ['E2', 'step2']
      _tmp = pd.concat([_term1, _term2], axis=1)\
               [['step','step2']].max(axis=1)
      _tmp.columns = ['step']
      _out['step'] = _tmp
      _term2.columns = ['E', 'step']
      _term3 = _term1 / _term2
      _out['E'] = _term3['E']
    elif isinstance(other, float) or isinstance(other, int):
      _out['E'] = self['E'] / float(other)
    return _out

  ###################
  # unit conversion #
  ###################
  # NOTE: energy MUST be read in as Hartree
  def ev(self):
    return self * 27.21138505

  def kcal(self):
    return self * 627.509469

class ScatterPlot(object):
  def __init__(self, QMData_pred, QMData_true, unit):
    # restructure data
    self._pred = pd.DataFrame(QMData_pred)
    self._pred.columns = ['Epred', 's']
    self._true = pd.DataFrame(QMData_true)
    self._true.columns = ['Etrue', 'step']
    self._step = 500

    # select only both valid and small steps
    self.scatter = pd.concat([self._pred, self._true], axis=1)
    self.scatter = self.scatter[np.isfinite(self.scatter['Epred'])]
    self.scatter = self.scatter[np.isfinite(self.scatter['Etrue'])]
    self.scatter = self.scatter[self.scatter.step < self._step]
    del self.scatter['s']

    # construct linear regression
    self._x = self.scatter['Epred'].values
    self._y = self.scatter['Etrue'].values
    if len(self._x) == 0:
      sys.exit("ERROR from analysis.py->ScatterPlot: " +\
               "Empty data array! file name matched? " +\
               "use QMData.extract('regexPatter')"
              )
    self.fit = np.polyfit(self._x, self._y, 1)
    self._y_fit = self.fit[0]*self._x + self.fit[1]

    # plot data
    self.data = np.array([self._x, self._y, self._y_fit]).T
    self._diff = abs(self._y - self._y_fit)
    self.MAE = sum(self._diff)/len(self._diff)
    self.RMSE = np.sqrt(sum(np.square(self._diff))/len(self._diff))
    self.xlabel = '$E_{pred}$ [' + unit + ']'
    self.ylabel = '$E_{true}$ [' + unit + ']'

  def plot(self):
    plot_data = self.data
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
    ax.set_xlabel(self.xlabel, fontsize=20)
    ax.set_ylabel(self.ylabel, fontsize=20)
    ax.tick_params(labelsize=15)

  def show(self):
    plt.show()

  def save(self, out_file):
    self.fig.savefig(out_file)

class DensityPlot(object):
  def __init__(self, path):
    self.out_dir = QMResults(path)
