import qctoolkit as qtk
import re, sys, os, copy, shutil
import numpy as np
from qctoolkit import utilities as ut
import qctoolkit.io_format.setting_pw as pw
import qctoolkit.io_format.pwinp as qin

def qmDir(inp, **kwargs):
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
    qtk.warning("io_format.run_folder.qmDir_inplace:"+\
                " output file "+ qps(inpdir+psinp)+\
                '.out exist, nothing to be done')
    new_run = False

  return qps(inpdir), qps(inpname), qps(psinp), new_run, kwargs
