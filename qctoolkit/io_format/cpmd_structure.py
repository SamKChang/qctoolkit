import qctoolkit as qtk
import re, sys, os, copy, shutil
import numpy as np
from qctoolkit import utilities as ut
import qctoolkit.io_format.setting_pw as pw
import qctoolkit.io_format.pwinp as qin

def qmDir_inplace(inp, **kwargs):
  qps = qtk.pathStrip
  _prefix = ''
  if 'prefix' in kwargs:
    _prefix = kwargs['prefix']
    del kwargs['prefix']
  _suffix = ''
  if 'suffix' in kwargs:
    _suffix = kwargs['suffix']
    del kwargs['suffix']
  try:
    root = re.match(re.compile('(.*/)[^\/]*'),inp).group(1)+'/'
  except:
    root = './'

  inproot = re.sub('\.inp', '',re.sub('.*/', '', inp))
  psinp = _prefix + inproot + _suffix
  inpdir = root
  inpname = inpdir + psinp + ".inp"
  new_run = True
  if os.path.exists(inpdir+psinp+'.out'):
    qtk.warning("io_format.cpmd.qmDir_inplace: output file "+\
                qps(inpdir+psinp)+\
                '.out exist, nothing to be done')
    new_run = False

  return qps(inpdir), qps(inpname), qps(psinp), new_run, kwargs

def qmDir(inp, **kwargs):
  """
  an root/inp folder contains all inp files
  inp files in root/inp/foo.inp
  will be copied to root/inp/foo/foo.inp
  scratch files will be generated and cleaned at root/inp/fc
  """

  qps = qtk.pathStrip
  _prefix = ''
  if 'prefix' in kwargs:
    _prefix = kwargs['prefix']
    del kwargs['prefix']
  _suffix = ''
  if 'suffix' in kwargs:
    _suffix = kwargs['suffix']
    del kwargs['suffix']
  _inplace = False
  if 'outdir' in kwargs:
    outdir = re.sub('\/$','', kwargs['outdir']) + '/'
    del kwargs['outdir']
  else:
    outdir = './'
  try:
    root = re.match(re.compile('(.*/)[^\/]*'),inp).group(1)\
           + outdir
  except:
    root = './' + outdir

  inproot = re.sub('\.inp', '',re.sub('.*/', '', inp))
  psinp = _prefix + inproot + _suffix
  inpdir = root + psinp
  inpname = inpdir + "/" + psinp + ".inp"
  new_run = True
  if not os.path.exists(inpdir):
    os.makedirs(inpdir)
    shutil.copyfile(inp, inpname) # copy inp file to folder
  elif _inplace:
    shutil.copyfile(inp, inpname)
  else:
    qtk.warning("io_format.cpmd.qmDir: folder '" + inpdir +\
               "' exists, nothing to be done")
    new_run = False
  # return path names for inpdir, inpname, psinp
  return qps(inpdir), qps(inpname), qps(psinp), new_run, kwargs
