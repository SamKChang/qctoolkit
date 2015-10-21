import qctoolkit as qtk
import re, sys, os, copy, shutil
import numpy as np
from qctoolkit import utilities as ut
import qctoolkit.io_format.setting_pw as pw
import qctoolkit.io_format.pwinp as qin

def qmDir(inp, **kwargs):
  """
  take input file as string, 
  create working folder,
  return corresponding names
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
  if 'inplace' in kwargs and kwargs['inplace']:
    _inplace = True
  if 'outdir' in kwargs:
    outdir = re.sub('\/$','', kwargs['outdir'])
    del kwargs['outdir']
  else:
    outdir = ''

  try:
    root = os.path.join(
             re.match(re.compile('(.*/)[^\/]*'),inp).group(1),
             outdir)
  except:
    root = os.path.join('.', outdir)

  inproot = re.sub('\.inp', '',re.sub('.*/', '', inp))
  psinp = _prefix + inproot + _suffix
  if _inplace:
    inpdir = root
  else:
    inpdir = os.path.join(root, psinp)
  inpname = os.path.join(inpdir,  psinp) + ".inp"
  new_run = True
  if _inplace:
    if os.path.exists(inpdir+psinp+'.out'):
      qtk.warning("io_format.run_folder.qmDir_inplace:"+\
                  " output file "+ qps(inpdir+psinp)+\
                  '.out exist, nothing to be done')
  else:
    if not os.path.exists(inpdir):
      os.makedirs(inpdir)
      shutil.copyfile(inp, inpname) # copy inp file to folder
  #    elif _inplace:
  ##      shutil.copyfile(inp, inpname)
    else:
      qtk.warning("io_format.run_folder.qmDir: folder '" + \
                  inpdir + "' exists, nothing to be done")
  #    new_run = False
  #  # return path names for inpdir, inpname, psinp
  return qps(inpdir), qps(inpname), qps(psinp), new_run, kwargs
