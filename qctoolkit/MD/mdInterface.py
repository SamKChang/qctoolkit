import qctoolkit as qtk
import mdcode.cpmd as cpmd
import os

def MDOut(out_folder, **kwargs):
  if 'program' not in kwargs:
    kwargs['program'] = qtk.setting.mdcode
  if kwargs['program'].lower() == 'cpmd':
    return cpmd.out(out_folder, **kwargs)
