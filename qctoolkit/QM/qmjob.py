import multiprocessing as mp
import qmout
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
  if 'save_restart' in kwargs:
    _save_restart = kwargs['save_restart']
  else:
    _save_restart = False

  ###########################################
  # SYSTEM CALL SUBPROCESS: Running mpi job #
  ###########################################
  # switch to working directory for each calculation
  # MUST switch pack by calling 'os.chdir(cwd)' at the end
  def compute(exestr, outpath, threads_per_job):
    outfile = open(outpath, "w")
    run = sp.Popen("%s %d %s %s"\
                   % (setting.mpistr, threads_per_job,
                      setting.ompstr, exestr), 
                   shell=True,
                   stdout=outfile)
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

    out = os.path.splitext(inp)[0] + '.out'
    exestr = "%s %s" % (exe, inp)
    compute(exestr, out, _threads)

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
      run.wait()
      log.close()
  
    if os.path.exists(out):
      qio_out = qmout.QMOut(out, program)
    else:
      qio_out = None

    return qio_out

  #######################
  # VASP IMPLEMENTATION #
  #######################
  elif program.lower() == 'vasp':
    if 'exe' in kwargs:
      exestr = kwargs['exe']
    else:
      exestr = setting.vasp_exe
    qmoutput = inp + '.out'
    compute(exestr, qmoutput, _threads)
    qio_out = qmout.QMOut('vasprun.xml', program)

    if not _save_restart:
      try:
        os.remove('WAVECAR')
      except:
        pass
    try:
      os.remove('POTCAR')
    except:
      pass

    if _delete:
      shutil.rmtree(inpdir)

    return qio_out

  # !!!!! TODO LIST !!!!! #
  #############################
  # Gaussian09 IMPLEMENTATION #
  #############################
  elif program.lower() == 'gaussian':
    ut.exit("ERROR! program '%s' not implemented" % program)
  ##################################
  # QuantumESPRESSO IMPLEMENTATION #
  ##################################
  elif program.lower() == 'pwscf':
    ut.exit("ERROR! program '%s' not implemented" % program)
  #########################
  # NwChem IMPLEMENTATION #
  #########################
  elif program.lower() == 'nwchem':
    ut.exit("ERROR! program '%s' not implemented" % program)
  #########################
  # GAMESS IMPLEMENTATION #
  #########################
  elif program.lower() == 'gamess':
    ut.exit("ERROR! program '%s' not implemented" % program)
  # and others... #

  else: 
    ut.exit("ERROR! program '%s' not recognized" % program)
