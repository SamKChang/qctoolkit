import qctoolkit as qtk
from xyzio import xyzOutput
import os

def MDOut(traj):
  stem, ext = os.path.splitext(traj)
  if ext == 'xyz':
    return xyzOutput(traj)
