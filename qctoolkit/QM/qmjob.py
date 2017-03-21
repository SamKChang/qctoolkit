import qctoolkit as qtk
import multiprocessing as mp
import qmInterface as qio
import glob, re, os, shutil
import fileinput
import subprocess as sp
import qctoolkit.utilities as ut
import numpy as np
import qctoolkit.setting as setting

# python interface to run QM code
# all code dependent part should be wrapped here
# inp is passed as relative path to the calling script path
def QMRun(inp, program=setting.qmcode, **kwargs):
  """
  interface to run qmcode with written inp files. It manages to start
  qmcode in **cwd**, write to output, collect result, 
  clean up scratch, tmp files according to setup. 
  However, each code require different implementation.

  input:
    inp(str): path to input file
    code(str): qmcode to run, default set to setting.qmcode

  kwargs (optional):
    threads=n(int): number of threads per job
    bigmem=Boolean: big memory request, implemented for CPMD and others

    CPMD:
      save_restart=Boolean
      scr=/path/to/scratch
  """

  if 'threads' in kwargs:
    _threads = kwargs['threads']
  else:
    _threads = 1
  if 'bigmem' in kwargs:
    _bigmem = kwargs['bigmem']
  else:
    if setting.memory > 16:
      _bigmem = True
    else:
      _bigmem = False
  if 'omp' in kwargs:
    omp = kwargs['omp']
  else:
    omp = 1
  if 'save_restart' in kwargs:
    _save_restart = kwargs['save_restart']
  else:
    _save_restart = False
  if 'save_density' in kwargs:
    _save_density = kwargs['save_density']
  else:
    _save_density = False

  if 'chdir' in kwargs and kwargs['chdir']:
    cwd = os.getcwd()
    os.chdir(inp)

  ###########################################
  # SYSTEM CALL SUBPROCESS: Running mpi job #
  ###########################################
  def compute(exestr, outpath, threads_per_job, **kwargs):
    """
    initiate a single MPI job, wait for it, and write to output
    """

    os.environ["OMP_NUM_THREADS"] = str(omp)

    outfile = open(outpath, "w")
    if threads_per_job > 1:
      mpi_cmd = "%s %d"% (setting.mpistr, threads_per_job)
      for mpi_flag in setting.mpi_flags:
        if mpi_flag == '--cpus-per-proc':
          flag = mpi_flag + ' ' + str(threads_per_job)
        else:
          flag = mpi_flag
        mpi_cmd = mpi_cmd + ' ' + flag
      cmd = mpi_cmd + ' ' + exestr
    else:
      cmd = exestr
    ut.progress('QMInp.run', 'running job with command: %s\n' % cmd)
    run = sp.Popen(cmd, shell=True, stdout=outfile)
    # wait each mpijob to finish before lauching another
    # otherwise all mpijobs will be launched simutaniously
    run.wait()
    outfile.close()
  ########## END OF SYSTEM CALL ##########

  #######################
  # CPMD IMPLEMENTATION #
  #######################
  if program.lower() == 'cpmd':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.cpmd_exe

    if 'scr' in kwargs and kwargs['scr']:
        scrdir = kwargs['scr']
        ut.delete(inpname, 'FILEPATH', 2)
        ut.insert(inpname, 'CPMD', ' FILEPATH\n  %s' % scrdir)

    inp_list = sorted(glob.glob('*.inp'))
    for job in inp_list:
      out = os.path.splitext(job)[0] + '.out'
      exestr = "%s %s" % (exe, job)
      compute(exestr, out, _threads)
      rst_list = sorted(glob.glob('RESTART.*'))
      if rst_list:
        rst_n = rst_list[-1]
        if os.path.exists('RESTART'):
          os.remove('RESTART')
        os.link(rst_n, 'RESTART')
    final_out = re.sub(r'_[0-9][0-9].out$', '.out', out)
    if final_out != out:
      os.copy(out, final_out)

    # clean up files
    files = sorted(glob.glob('*'))
    tmp = filter(\
      lambda x: '.out' not in x \
                and '.inp' not in x\
                and '.psp' not in x\
                and '.xyz' not in x\
                and 'KPTS_GENERATION' not in x\
                and 'RESTART' not in x\
                and 'DENSITY' not in x\
                and 'SPINDEN' not in x, files
    )
    for f in tmp: os.remove(f)
    if not _save_restart:
      rst_list = sorted(glob.glob("RESTART*"))
      for rfile in rst_list:
        os.remove(rfile)

    densities = sorted(glob.glob('*DEN*'))
    for i in range(len(densities)):
      exe = setting.cpmd_cpmd2cube_exe
      log_name = densities[i] + '_%02d.log' % i
      log = open(log_name, 'w')
      run = sp.Popen("%s -fullmesh %s" % (exe, densities[i]), 
               shell=True,
               stdout=log)
      run.wait()
      log.close()
  
    if os.path.exists(out):
      qio_out = qio.QMOut(out, program='cpmd')
    else:
      qio_out = None

  #######################
  # VASP IMPLEMENTATION #
  #######################
  elif program.lower() == 'vasp':
    if 'exe' in kwargs:
      exestr = kwargs['exe']
    else:
      exestr = setting.vasp_exe
    qmoutput = inp + '.out'
    qmlog = inp + '.log'
    compute(exestr, qmoutput, _threads)
    qio_out = qio.QMOut('vasprun.xml', program='vasp')

    if not _save_restart:
      try:
        os.remove('WAVECAR')
      except:
        pass
    try:
      os.remove('POTCAR')
    except:
      pass

    os.rename(qmoutput, qmlog)
    shutil.copyfile('vasprun.xml', qmoutput)

  #########################
  # NWChem IMPLEMENTATION #
  #########################
  elif program.lower() == 'nwchem':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.nwchem_exe
    exestr = "%s %s" % (exe, inp)
    qmoutput = os.path.splitext(inp)[0] + '.out'
    compute(exestr, qmoutput, _threads)
    qio_out = qio.QMOut(qmoutput, program='nwchem')

    files = sorted(glob.glob('*.*'))
    tmp = filter(\
      lambda x: '.out' not in x \
                and '.inp' not in x\
                and '.cube' not in x\
                and '.movecs' not in x, files
    )
    for f in tmp: os.remove(f)
    movecs = sorted(glob.glob('*.movecs'))
    for f in movecs:
      exe = setting.nwchem_mov2asc_exe
      nb = qio_out.n_basis
      out = re.sub('\.movecs','.modat',f)
      exestr = "%s %d %s %s" % (exe, nb, f, out)
      run = sp.Popen(exestr, shell=True)
      run.wait()
      qio_out.getMO(out)
      if not _save_restart:
        os.remove(f)

  #########################
  # BigDFT IMPLEMENTATION #
  #########################
  elif program.lower() == 'bigdft':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.bigdft_exe
    inp = os.path.splitext(inp)[0]
    exestr = "%s %s" % (exe, inp)
    qmoutput = inp + '.out'
    compute(exestr, qmoutput, _threads)
    qio_out = qio.QMOut(qmoutput, program='bigdft')

  #########################
  # Abinit IMPLEMENTATION #
  #########################
  elif program.lower() == 'abinit':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.abinit_exe
    inp = os.path.splitext(inp)[0]
    exestr = "%s < %s" % (exe, inp + '.files')
    qmlog = inp + '.log'
    qmoutput = inp + '.out'
    compute(exestr, qmlog, _threads)

    # clean up files
    files = sorted(glob.glob('*'))
    tmp = filter(\
      lambda x: '.out' not in x \
                and '.log' not in x\
                and '.inp' not in x\
                and '.files' not in x\
                and 'psp' not in x\
                and '_EIG' not in x\
                and '_WFK' not in x\
                and '_DOS' not in x\
                and '_DEN' not in x, files
    )
    for f in tmp: os.remove(f)
    if not _save_restart:
      for rst in sorted(glob.glob('*_WFK')):
        os.remove(rst)
    if not _save_density:
      for den in sorted(glob.glob('*_DEN')):
        os.remove(den)

    densities = sorted(glob.glob('*_DEN'))
    if len(densities) > 0:
      i = len(densities) - 1
      exe = setting.abinit_cut3d_exe
      #den_inp = densities[i] + '.cov'
      den_inp = densities[-1] + '.cov'
      cube_file = inp + '.cube'
      den_inp_file = open(den_inp, 'wb')
      den_inp_file.write("%s\n1\n14\n%s" % (densities[i], cube_file))
      den_inp_file.close()
      log_name = densities[i] + '.log'
      log = open(log_name, 'w')
      run = sp.Popen("%s < %s" % (exe, den_inp), 
               shell=True,
               stdout=log)
      run.wait()
      log.close()

    qio_out = qtk.QMOut(qmoutput, program='abinit')

  #############################
  # Gaussian09 IMPLEMENTATION #
  #############################
  elif program.lower() == 'gaussian':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.gaussian_exe

    exestr = "%s %s" % (exe, inp)
    qmoutput = os.path.splitext(inp)[0] + '.out'
    qmlog = os.path.splitext(inp)[0] + '.log'
    compute(exestr, qmoutput, _threads)
    os.rename(qmlog, qmoutput)

    chks = sorted(glob.glob('*.chk'))
    for chk in chks:
      exe = setting.gaussian_formchk_exe
      exestr = "%s %s" % (exe, chk)
      run = sp.Popen(exestr, shell=True)
      run.wait()
    if _save_density:
      fchk = sorted(glob.glob("*.fchk"))[-1]
      q = qtk.CUBE(fchk)

    qio_out = qio.QMOut(qmoutput, program='gaussian')

  ##################################
  # QuantumESPRESSO IMPLEMENTATION #
  ##################################
  elif program.lower() == 'espresso':
    if 'exe' in kwargs:
      exe = kwargs['exe']
    else:
      exe = setting.espresso_exe

    inp_list = sorted(glob.glob('*.inp'))
    for job in inp_list:
      out = os.path.splitext(job)[0] + '.out'
      exestr = "%s < %s" % (exe, job)
      compute(exestr, out, _threads)
    qio_out = qio.QMOut(out, program='espresso')
    if not _save_restart:
      rst_list = sorted(glob.glob("*.wfc*"))
      rst_list.extend(sorted(glob.glob("*.restart_*")))
    else:
      rst_list = []
    for r in rst_list:
      os.remove(r)
      
  #########################
  # GAMESS IMPLEMENTATION #
  #########################
  elif program.lower() == 'gamess':
    ut.exit("ERROR! program '%s' not implemented" % program)
  # and others... #

  else: 
    ut.exit("ERROR! program '%s' not recognized" % program)

  if 'cwd' in locals():
    os.chdir(cwd)

  qio_out.path = inp
  return qio_out
