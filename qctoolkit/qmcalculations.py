import multiprocessing as mp
import qmio as qio
import glob, re, os, shutil
import fileinput
import subprocess as sp
import utilities as ut
import numpy as np
import setting

# python interface to run QM code
# all code dependent part should be wrapped here
# inp is passed as relative path to the calling script path
def QMRun(inp, program, **kwargs):
  if 'threads' in kwargs:
    _threads = kwargs['threads']
  else:
    _threads = 1
  _delete = False
  if 'QMReturn' in kwargs:
    # return value of QMRun
    # return = [list of properties to return]
    _return = kwargs['QMReturn']
  else:
    _return = False
  if 'cleanup' in kwargs and kwargs['cleanup']:
    _delete = True

  if 'bigmem' in kwargs:
    _bigmem = kwargs['bigmem']
  else:
    if mp.cpu_count() > 20:
      _bigmem = True
    else:
      _bigmem = False
  if 'prefix' in kwargs:
    _prefix = kwargs['prefix']
  else:
    _prefix = ''
  if 'suffix' in kwargs:
    _suffix = kwargs['suffix']
  else:
    _suffix = ''

  ###########################################
  # SYSTEM CALL SUBPROCESS: Running mpi job #
  ###########################################
  # switch to working directory for each calculation
  # MUST switch pack by calling 'os.chdir(cwd)' at the end
  def compute(exestr, outpath, threads_per_job):
    outfile = open(outpath, "w")
    run = sp.Popen("mpirun -np %d %s"\
                   % (threads_per_job, exestr), 
                   shell=True,
                   stdout=outfile)
    # wait each mpijob to finish before lauching another
    # otherwise all mpijobs will be launched simutaniously
    run.communicate() 
    outfile.close()
  ########## END OF SYSTEM CALL ##########

  #######################
  # CPMD IMPLEMENTATION #
  #######################
  # an root/inp folder contains all inp files
  # inp files in root/inp/foo.inp
  # will be copied to root/inp/foo/foo.inp
  # scratch files will be generated and cleaned at root/inp/fc
  if program == 'cpmd':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.cpmd_exe
    cwd = os.getcwd()
    if 'save_restart' in kwargs:
      _save_restart = True
    else:
      _save_restart = False

    # I/O setup for file names, folder etc
    # out dir is set to be relative to inp file
    # default: same as inp file
    if 'outdir' in kwargs:
      outdir = re.sub('\/$','', kwargs['outdir']) + '/'
    else:
      outdir = './'
    try:
      root = re.match(re.compile('(.*/)[^\/]*'),inp).group(1)\
             + outdir
    except:
      root = './' + outdir
    if 'alchemRef' in kwargs:
      _alchemRef = True
      inproot = "rst_" + re.sub('\.inp', '',re.sub('.*/', '', inp))
    else:
      _alchemRef = False
      inproot = re.sub('\.inp', '',re.sub('.*/', '', inp))
    inpdir = root + _prefix + inproot + _suffix
    inpname = inpdir + "/" + inproot + ".inp"
    print inpdir, inpname

    # create new dir for CPMD calculation
    if not os.path.exists(inpdir):
      os.makedirs(inpdir)
      shutil.copyfile(inp,inpname) # copy inp file to folder
    else:
      return np.nan
      ut.exit("folder '" + inpdir +\
              "' exist, nothing to be done")

    if 'scr' in kwargs and kwargs['scr']:
      _scr = True
      scrdir = kwargs['scr']
    else:
      _scr = False
      scrdir = inpdir

    # RESTART file for alchemy
    if 'alchemScan' in kwargs and kwargs['alchemScan']:

      if 'alchemRefPrefix' in kwargs:
        alchemPrefix = kwargs['alchemRefPrefix']
      else:
        alchemPrefix = 'rst_'

      if not 'alchemRefPath' in kwargs:
        ut.exit("alchemical reference not set")
      elif not os.path.exists(kwargs['alchemRefPath']):
        ut.exit("alchemical reference '%s' not found"\
                % kwargs['alchemRefPath'])
      else:
        _alchemScan = True
        ref = kwargs['alchemRefPath']
        refroot = alchemPrefix + re.sub('\.inp', '',
                                  re.sub('.*/', '', ref))

        # for flexibility

        # for alchemy QMJobs
        #refpath = re.sub(inproot, refroot, scrdir)

        # !!!!BUG!!!! alchemy QMJobs need to be modified!!!
        # for general interface
        refpath = ref
        rst_src = refpath + "/RESTART.1"
        rst_trg = inpdir + "/RESTART"
        os.link(rst_src, rst_trg)

    else:
      _alchemScan = False

    ####################################################
    # PROCESS INP FILES: according to calculation mode #
    ####################################################
    # supress writing RESTART
    if ('debug' in kwargs and kwargs['debug']) or _return:
      ut.insert(inpname, 
                "CPMD", 
                " BENCHMARK\n  1 0 0 0 0 0 0 0 0 0")

    # alchemical prediction
    if _alchemScan:
      ut.delete(inpname,'MAXITER', 2)
      ut.insert(inpname,'MIRROR', ' MAXITER\n  1')
      ut.delete(inpname, 'RESTART', 1)
      ut.delete(inpname, 'INITIALIZE WAVEFUNCTION', 1)
      ut.insert(inpname,'CPMD', ' RESTART WAVEFUNCTION')

    if _alchemRef:
      ut.delete(inpname, 'BENCHMARK', 2)
      ut.delete(inpname, 'RESTART', 1)

    # big memory setup
    #if 'bigmem' in kwargs:
    if _bigmem:
      ut.delete(inpname, 'MEMORY BIG', 1)
      ut.insert(inpname, 'MIRROR', ' MEMORY BIG')
    else:
      ut.delete(inpname, 'MEMORY BIG', 1)

    # scratch path
    if _scr:
      ut.delete(inpname, 'FILEPATH', 2)
      ut.insert(inpname, 'CPMD', ' FILEPATH\n  %s' % scrdir)

    # !!!! TODO !!!! #
    # scanning alchemical path
    # second order alchemy
    # KSEg
    # geopt
    # metal surface
    # different optimizer?
    ########## END OF PROCESSING INP FILES ##########

    os.chdir(inpdir)
    exestr = "%s %s" % (exe, inproot + ".inp")
    qmoutput = inproot + ".out"
    compute(exestr, qmoutput, _threads)

    # clean up files
    try:
      os.remove('LATEST')
    except OSError:
      pass
    try:
      os.remove('GEOMETRY.scale')
    except OSError:
      pass
    try:
      os.remove('GEOMETRY')
    except OSError:
      pass
    try:
      os.remove('KPTS_GENERATION')
    except OSError:
      pass
    try:
      os.remove('RESTART')
    except OSError:
      pass
    if not _save_restart:
      rst_list = glob.glob("RESTART*")
      for rfile in rst_list:
        os.remove(rfile)

    if os.path.exists('DENSITY'):
      exe = setting.cpmd_cpmd2cube
      log = open('DENSITY.log', 'w')
      run = sp.Popen("%s -fullmesh DENSITY" % exe, 
               shell=True,
               stdout=log)
      run.communicate()
      log.close()
      

    qio_out = qio.QMOut(qmoutput, program)
    os.chdir(cwd)
    
    if _delete:
      shutil.rmtree(inpdir)

  # !!!!! TODO LIST !!!!! #
  #############################
  # Gaussian09 IMPLEMENTATION #
  #############################
  elif program == 'gaussian':
    ut.exit("ERROR! program '%s' not implemented" % program)
  ##################################
  # QuantumESPRESSO IMPLEMENTATION #
  ##################################
  elif program == 'pwscf':
    ut.exit("ERROR! program '%s' not implemented" % program)
  #########################
  # NwChem IMPLEMENTATION #
  #########################
  elif program == 'nwchem':
    ut.exit("ERROR! program '%s' not implemented" % program)
  #########################
  # GAMESS IMPLEMENTATION #
  #########################
  elif program == 'gamess':
    ut.exit("ERROR! program '%s' not implemented" % program)
  # and others... #

  else: 
    ut.exit("ERROR! program '%s' not recognized" % program)

  return qio_out


# collection of QM input files, run all in parallel
# QM code independent implementation
# target folder should contain many inp files for QM scan
class QMJobs(object):
  def __init__(self, path, pattern, program, **kwargs):

    self._threadspj = 1
    self._path = re.sub(re.compile('/$'), '', path)
    self._inps = sorted(glob.glob(self._path + "/" + pattern))
    self._program = program
    self._alchemScan = False

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

    if 'scr' in kwargs:
      self._scr = kwargs['scr']
    else:
      self._scr = False

    if 'debug' in kwargs:
      self._debug = kwargs['debug']
    else:
      self._debug = False
   
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
  # one more layer to organize different running mode
  # calling QMRun()
  def single_run(self, inp, **kwargs):
    if 'threads' in kwargs:
      _job_threads = kwargs['threads']
    else:
      _job_threads = 1

    if self._alchemScan:
      if 'alchemRef' in kwargs and kwargs['alchemRef']:
        QMRun(inp, self._program, 
              debug = self._debug,
              scr=self._scr,
              threads = _job_threads,
              save_restart=True,
              outdir='../',
              alchemRef=True)
      else:
        QMRun(inp, self._program, 
              debug = self._debug,
              outdir='../',
              threads = _job_threads,
              scr=self._scr,
              alchemRefPath=self._ref,
              alchemScan=True)
    else:
      QMRun(inp, self._program, 
            debug = self._debug,
            outdir='../', 
            threads = _job_threads,
            scr=self._scr)

  # BUGS!!!! JoinableQueue might be necessary!      
  #########################
  # MAIN PARALLEL ROUTINE #
  #########################
  # reference scan of true energy
  # QM code independent
  def run(self, **kwargs):
    # default, 1 core per job, all cores used
    if 'threads_per_job' in kwargs:
      self._threadspj = kwargs['threads_per_job']
    else:
      self._threadspj = 1
    if 'threads' in kwargs:
      self._threads = kwargs['threads']
    else:
      self._threads = mp.cpu_count()

    # subroutine for running queued jobs
    # it runs until the signal 'DONE' is reached
    # calling single_run
    def run_jobs(inp_queue):
      inp = inp_queue.get()
      while (inp != 'DONE'):
        # run single job
        self.single_run(inp,
                        threads=self._threadspj)
        inp = inp_queue.get()
      if (inp == 'DONE'): 
        # put 'DONE' back to the queue
        # it's necessary for other run_jobs process
        # to get 'DONE' signal
        inp_queue.put('DONE')
    # end of runing queued jobs subroutine #

    if self._alchemScan:
      ut.report("Alchemy scan", "running reference", 
                self._ref, "with", self._threads, "cores")
      self.single_run(self._ref, 
                      alchemRef=True,
                      threads=self._threads)
      ut.report("Alchemy scan", "reference calculation done")

    #jobs = []
    inp_queue = mp.Queue()
    for inp in self._inps: # construct inp queue
      inp_queue.put(inp)
    inp_queue.put('DONE')  # attach 'DONE' to the end of queue
    for core in range(self._threads): # start run_job on each core
      p = mp.Process(target=run_jobs,\
          args=(inp_queue, ))
      p.start()
    #  jobs.append(p)
    #for process in jobs: # wait for all jobs to finish
    #  process.join()    
  ####### END OF MAIN PARALLEL ROUTIN #######
