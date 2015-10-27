import qctoolkit.setting as setting
import multiprocessing as mp
import operator
from compiler.ast import flatten

def parallelize(target_function, 
                input_list, 
                **kwargs):
 
  if 'threads' in kwargs:
    threads = kwargs['threads']
  else:
    threads = setting.cpu_count
  if 'block_size' in kwargs:
    block_size = kwargs['block_size']
  else:
    block_size = len(input_list)/(threads*3)

  #############################################
  # runing target function of a single thread #
  #############################################
  def run_jobs(q_in, q_out):
    for inps in iter(q_in.get, None):
      ind = inps[-1]    # index of job
      inps = inps[:-1]  # actual input sequence
      out = []
      for args in inps:
        if type(args[-1]) == dict: # check known args input
          kwargs = args[-1]
          args = args[:-1]  
          out.append(target_function(*args, **kwargs))
        else:
          out.append(target_function(*args))
      if out != None:
        q_out.put([out, ind]) # output result with index
      q_in.task_done()
    q_in.task_done() # task done for 'None' if q_in finished
  ###### end of single thread definition ######

  # devide input_list into chunks according to block_size
  def chunks(_list, _size):
    for i in range(0, len(_list), _size):
      yield _list[i:i+_size]
  input_block = list(chunks(input_list, block_size))

  # setup empty queue
  output_stack = []
  output = []
  qinp = mp.JoinableQueue()
  qout = mp.Queue()

  # start process with empty queue
  for thread in range(threads):
    p =  mp.Process(target=run_jobs, args=(qinp, qout))
    p.daemon = True # necessary for terminating finished thread
    p.start()

  # put I/O data into queue for parallel processing
  index = range(len(input_block))
  for ind, inps in zip(index, input_block):
    inps.append(ind) # append inp index
    qinp.put(inps)   # put inp to input queue
  qinp.join()       # wait for jobs to finish

  # 'while not queue.empty' is NOT reliable
  if not qout.empty():
    for i in range(len(input_block)):
      output_stack.append(qout.get())

  if len(output_stack)>0:
    # sort/restructure output according to input order
    output_stack = sorted(output_stack, key=operator.itemgetter(1))
    # loop though all input for corresponding output
    for data_out in output_stack: 
      # if output is list of class, in-line iteration doesn't work
      output.append(data_out[0])
    return flatten(output)
