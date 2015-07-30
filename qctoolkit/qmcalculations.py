import multiprocessing as mp
import glob, re, os, shutil
import fileinput
import subprocess as sp
import utilities as ut

# python interface to run QM code
# all code dependent part should be wrapped here
def QMRun(inp, program, **kwargs):
  if 'threads_per_job' in kwargs:
    _threads = kwargs['threads_per_job']
  else:
    _threads = 1

  #######################
  # CPMD IMPLEMENTATION #
  #######################
  if program == 'cpmd':
    exe = 'cpmd.x'
    if 'save_restart' in kwargs:
      _save_restart = True
    else:
      _save_restart = False

    # I/O setup 
    inproot = re.sub('.*/', '', inp)
    inpdir = re.sub('/'+inproot, '', inp)
    cwd = os.getcwd()

    if 'alchemScan' in kwargs:
      # prepare RESTART file
      ref = kwargs['alchemRef']
      ref_path = re.sub('inp\/', 'rst_', ref)
      ref_path = re.sub('\.inp','/RESTART.1', ref_path)
      rst_path = inpdir + "/RESTART"
      os.link(ref_path, rst_path)
      # modify inp file
      ut.delete(inp,'MAXITER', 2)
      ut.insert(inp,'MIRROR', ' MAXITER\n  1')
      ut.delete(inp, 'RESTART', 1)
      ut.insert(inp,'CPMD', ' RESTART WAVEFUNCTION')

    if 'bigmem' in kwargs:
      _bigmem = kwargs['bigmem']
    else:
      _bigmem = False
    if _bigmem:
      ut.delete(inp, 'MEMORY BIG', 1)
      ut.insert(inp, 'MIRROR', ' MEMORY BIG')

    os.chdir(inpdir)
    out = re.sub('\.inp','.out', inproot)
    outfile = open(out, "w")
    run = sp.Popen("mpirun -np %d %s %s"\
                   % (_threads, exe, inproot), 
                   shell=True,
                   stdout=outfile)
    run.wait()

    # clean up files
    try:
      os.remove('LATEST')
      os.remove('GEOMETRY')
      os.remove('KPTS_GENERATION')
      os.remove('RESTART')
    except OSError:
      pass
    if not _save_restart:
      try:
        os.remove('RESTART.1')
      except OSError:
        pass

    os.chdir(cwd)


# collection of QM input files, run all in parallel
# QM code independent implementation
class QMJobs(object):
  def __init__(self, path, pattern, program, **kwargs):

    self._threadspj = 1
    self._path = re.sub(re.compile('/$'), '', path)
    self._inps = sorted(glob.glob(self._path + "/" + pattern))
    self._program = program

    if 'bigmem' in kwargs:
      self._bigmem = kwargs['bigmem']
    else:
      if mp.cpu_count() > 20:
        self._bigmem = True
      else:
        self._bigmem = False

    if 'outdir' in kwargs:
      self._outdir = kwargs['outdir']
    else:
      self._outdir = re.sub('$', '/../', self._path)
   
  # set alchemy scanning mode and define alchemical reference inp
  def setAlchemScan(self, **kwargs):
    self._alchemScan = True
    if 'ref' in kwargs:
      self._ref = self._path + "/" + kwargs['ref']
    else:
      self._ref = self._inps[0]

    if self._ref in self._inps:
      del self._inps[self._inps.index(self._ref)]
    else:
      ut.exit("reference file '" + self._ref + "' not exist")

  # run single job, called by parallel runner
  # calling QMRun()
  def single_run(self, inp, **kwargs):
    if 'threads' in kwargs:
      _job_threads = kwargs['threads']
    else:
      _job_threads = 1

    inproot = re.match(re.compile('.*/(.*)\..*'),inp).group(1)
    if 'alchem_ref' in kwargs:
      inpname = self._outdir\
                + "%s/%s.inp" % ("rst_"+inproot, inproot)
      inpdir = re.sub('$', "rst_"+inproot, self._outdir)
    else:
      inpname = self._outdir + "%s/%s.inp" % (inproot, inproot)
      inpdir = re.sub('$', inproot, self._outdir)

    if not os.path.exists(inpdir):
      os.makedirs(inpdir) # create working folder
      shutil.copyfile(inp,inpname) # copy inp file to folder
      if self._alchemScan:
        if 'alchem_ref' in kwargs:
          QMRun(inpname, self._program,
                threads_per_job=_job_threads,
                save_restart=True,
                bigmem=self._bigmem,
               )
        else:
          QMRun(inpname, self._program,
                threads_per_job=_job_threads,
                alchemScan=True,
                bigmem=self._bigmem,
                alchemRef=self._ref
               )
      else:
        QMRun(inpname, self._program,
              bigmem=self._bigmem,
              threads_per_job=_job_threads)

    # skip if working folder is already exist
    else:
      ut.warning("folder '" + inpdir +\
                 "' exist, nothing to be done")
      
  ################
  # MAIN ROUTINE #
  ################
  # reference scan of true energy
  # calling self.single_run()
  def run(self, **kwargs):
    if 'threads_per_job' in kwargs:
      self._threadspj = kwargs['threads_per_job']
    else:
      self._threadspj = 1
    if 'threads' in kwargs:
      self._threads = kwargs['threads']
    else:
      self._threads = mp.cpu_count()

    # run jobs parallel from queue
    def run_jobs(inp_queue):
      inp = inp_queue.get()
      while (inp != 'DONE'):
        # run single job
        self.single_run(inp,
                        threads=self._threadspj)
        inp = inp_queue.get()
      if (inp == 'DONE'):
        inp_queue.put('DONE')

    if self._alchemScan:
      ut.report("Alchemy scan: ", "running reference " + self._ref\
                + " with " + str(self._threads) + " cores")
      self.single_run(self._ref, 
                      alchem_ref=True,
                      threads=self._threads)
      ut.report("Alchemy scan: ", "reference calculation done")

    jobs = []
    inp_queue = mp.Queue()
    for inp in self._inps:
      inp_queue.put(inp)
    inp_queue.put('DONE')
    for core in range(self._threads):
      p = mp.Process(target=run_jobs,\
          args=(inp_queue, ))
      p.start()
      jobs.append(p)
    for process in jobs:
      process.join()    
