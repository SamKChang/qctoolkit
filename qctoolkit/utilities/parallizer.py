import qctoolkit as qtk
import qctoolkit.setting as setting
import multiprocessing as mp
import operator
from compiler.ast import flatten
import numpy as np
import sys, os
import shutil
import subprocess as sp
import pickle
import dill

# devide input_list into chunks according to block_size
def chunks(_list, _size):
  for i in range(0, len(_list), _size):
    yield _list[i:i+_size]

def qmWriteAll(inp_list, root, overwrite=False, 
  compress=False, file_per_folder=0):
  if os.path.exists(root):
    if overwrite:
      qtk.warning("overwrite existing folder %s" % root)
      shutil.rmtree(root)
      os.makedirs(root)
    else:
      qtk.warning("%s exists, " % root +\
        "joining calculations with other threads")
  else:
    os.makedirs(root)
  cwd = os.getcwd()
  os.chdir(root)
  if file_per_folder:
    inp_grp = list(chunks(inp_list, file_per_folder))
    n_digits = len(str(len(inp_grp)))
    for i, inps in enumerate(inp_grp):
      qmWriteAll(inps, 'jobs%s' % (str(i).zfill(n_digits)))
  else:
    for inp in inp_list:
      inp.write(inp.molecule.name)
  os.chdir(cwd)
  if compress:
    cmd = 'tar -zcf %s %s' % (root + '.tar.gz', root)
    run = sp.Popen(cmd, shell=True, stdin=sp.PIPE)
    run.stdin.flush()
    run.communicate()
    run.wait()
    qtk.report("qmWriteAll", "compression completed")

def qmRunAll(inp_list, root=None, **kwargs):
  if 'block_size' not in kwargs:
    kwargs['block_size'] = 1
  job = []
  for inp in inp_list:
    job.append([inp, inp.molecule.name])
  inp = inp_list[0]
  if inp.setting['threads'] != 1 and 'threads' not in kwargs:
    kwargs['threads'] = setting.cpu_count / inp.setting['threads']
  if root is None:
    qtk.parallelize(qtk.qmRunJob, job, **kwargs)
  else:
    if os.path.exists(root):
      if 'overwrite' in kwargs and kwargs['overwrite']:
        qtk.warning("overwrite existing folder %s" % root)
        shutil.rmtree(root)
        os.makedirs(root)
      else:
        qtk.warning("%s exists, " % root +\
          "joining calculations with other threads")
    else:
      os.makedirs(root)
    cwd = os.getcwd()
    os.chdir(root)
    qtk.parallelize(qtk.qmRunJob, job, **kwargs)
    os.chdir(cwd)

def qmRunJob(inp, name):
  qtk.report("qmRunJob", "runing qmjob:'%s'" % inp,
             'with name:', name)
  return inp.run(name)

def parallelize(target_function, 
                input_list, 
                n_output = 1,
                threads = None,
                **kwargs):
  """
  target_function is implemented in a general way
    supposely any function would work
    But it could break down if target_function assumes some convoluted data structure
    
  input_list is a list of list. 
    Each input entry should be wrapped properly as a list 
    **kwargs can be passed py passing dictionary
    
  Example:
    # a toy target function
    def f(a, b, **kwargs):
      if 'factor' in kwargs:
        factor = kwargs['factor']
      else:
        factor = 1
      return a + b*factor
      
    input_list = [[i,j,{'factor':3}] for i in range(10) for j in range(10)]
    
    out_list = parallelize(f, input_list, block_size=10)
  """

  if threads is None:
    threads = setting.cpu_count
  if 'block_size' in kwargs:
    block_size = kwargs['block_size']
  else:
    if len(input_list) > threads*3:
      block_size = len(input_list)/(threads*3)
    else:
      block_size = 1

  for inp in input_list:
    if not dill.pickles(inp):
      qtk.exit("pickling failed, subprocess only works if pickle is possible")

  #############################################
  # runing target function of a single thread #
  #############################################
  def run_jobs(q_in, q_out=None):
    for inps in iter(q_in.get, None):
      ind = inps[-1]    # index of job
      inps = inps[:-1]  # actual input sequence
      out = []
      try:
        for args in inps:
          if type(args[-1]) == dict: # check known args input
            kwargs = args[-1]
            args = args[:-1]  
            out.append(target_function(*args, **kwargs))
          else:
            out.append(target_function(*args))
        if q_out:
          q_out.put([out, ind]) # output result with index
      except: 
        qtk.warning('job failed!')
        if q_out:
          q_out.put([np.nan, ind])
  ###### end of single thread definition ######

  input_block = list(chunks(input_list, block_size))

  # setup empty queue
  output_stack = []
  output = []
  qinp = mp.Queue()
  qout = mp.Queue()

  # start process with empty queue
  jobs = []
  for thread in range(threads):
    if n_output > 0:
      p =  mp.Process(target=run_jobs, args=(qinp, qout))
    else:
      p =  mp.Process(target=run_jobs, args=(qinp,))
    p.start()
    jobs.append(p)

  # put I/O data into queue for parallel processing
  index = range(len(input_block))
  for ind, inps in zip(index, input_block):
    inps.append(ind) # append inp index
    qinp.put(inps)   # put inp to input queue

  for thread in jobs:
    qinp.put(None)

  # while not queue.empty' is NOT reliable
  for i in range(len(input_block)):
    # collect output from each subprocess
    try:
      if n_output > 0:
        output_stack.append(qout.get())
    # check keyboard interrupt and terminate subprocess
    except KeyboardInterrupt:
      qtk.warning('jobs terminated by keyboard interrupt')
      for p in jobs:
        p.terminate()
      try:
        sys.exit(0)
      except SystemExit:
        os._exit(0)

  for thread in jobs:
    thread.join()

  # clean up queues
  while not qinp.empty():
    qinp.get()
  if n_output > 0:
    while not qout.empty():
      qout.get()

  if len(output_stack)>0:
    # sort/restructure output according to input order
    output_stack = sorted(output_stack, key=operator.itemgetter(1))
    # loop though all input for corresponding output
    for data_out in output_stack: 
      # if output is list of class, in-line iteration doesn't work
      output.append(data_out[0])
  out = list(qtk.flatten(output))
  if n_output == 1:
    return out
  elif n_output > 1:
    out_list = []
    for i in range(n_output):
      out_list.append(out[i::n_output])
    return tuple(out_list)
  else:
    return 0
