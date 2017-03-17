from molecule import *
from QM.qmInterface import *
from QM.qmjob import *
from QM.qmresult import *
from QM.ofdft.libxc_dict import xc_dict
from alchemy.aljob import *
from alchemy.alpath import *
from analysis import CUBE
from analysis import PCA
from analysis import QMList
from QM.pseudo.pseudo import *
from utilities import *
from setting import *
from ccs.ccs import CCS
from QM.general_io import GenericQMInput as QMInput
from QM.general_io import GenericQMOutput as QMOutput
from data.elements.element_list import ELEMENTS as element
import MD
import ML
import ccs
import QM
import optimization
import alchemy
import setting
import DB

import os
import copy_reg
import copy
import types
import pickle

# check for qtk setting
missing_files = []
paths = os.environ["PATH"].split(":")
code_pattern = re.compile('cpmd|bigdft|vasp|nwchem|espresso')
exe_pattern = re.compile('.*exe')
url_pattern = re.compile('.*url')
for dep in dir(setting):
  if code_pattern.match(dep) and not url_pattern.match(dep):
    not_found = True
    file_str = getattr(setting, dep)
    if exe_pattern.match(dep):
       itr = 0
       while not_found and itr < len(paths):
         path = paths[itr]
         test_path = os.path.join(path, file_str)
         itr = itr +1
         if os.access(test_path, os.X_OK):
           not_found = False
    else:
      if os.access(file_str, os.F_OK):
        not_found = False
    if not_found:
      missing_files.append(file_str)

if missing_files:
  for missing_file in missing_files:
    qtk.warning("missing file: %s" % missing_file)
  qtk.warning(
               'please modify /path/to/qctoolkit/setting.py ' +\
               'and recompile.'
             )

# Steven Bethard's fix for instance method pickling
def _pickle_method(method):
  func_name = method.im_func.__name__
  obj = method.im_self
  cls = method.im_class
  return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
  for cls in cls.mro():
    try:
      func = cls.__dict__[func_name]
    except KeyError:
      pass
    else:
      break
  return func.__get__(obj, cls)

copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)
