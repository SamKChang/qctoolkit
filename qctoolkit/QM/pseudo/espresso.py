import qctoolkit as qtk
import re, os
import urllib2
import numpy as np
from qctoolkit.QM.general_io import InpContent
from cpmd import write as cpmd_write
import subprocess as sp

xc_dict = {
            'pbe': 1134,
            'lda': 900,
            'blyp': 1312,
            'factor': 0.6666666667,
            'bp': 1111,
          }

def write(self, cpmd_name, espresso_name):
  if not os.path.exists(cpmd_name):
    qtk.report("PP", "writing cpmd PP file")
    cpmd_write(self, cpmd_name)
    cpmd_exists = False
  else:
    cpmd_exists = True
    qtk.prompt('cpmd pp path:%s exist' % cpmd_name)
  if (cpmd_name == espresso_name and not cpmd_exists)\
  or not os.path.exists(espresso_name):
    qtk.report("PP", 'start converting Goedecker PP')
    conv_pp = sp.Popen("%s %s" % \
      (qtk.setting.espresso_cpmd2upf_exe, cpmd_name),
      shell=True)
    conv_pp.wait()
    if conv_pp.returncode != 0:
      qtk.warning('conversion failed...')
    else:
      os.rename(cpmd_name + '.UPF', espresso_name)
  else:
    qtk.prompt('espresso pp path:%s exist' % espresso_name)
