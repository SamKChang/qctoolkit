import qctoolkit as qtk
import numpy as np

def getCelldm(self):
  if 'celldm' not in self.setting:
    if not self.molecule.celldm:
      size = self.molecule.getSize()
      if 'margin' not in self.setting:
        m = qtk.setting.box_margin
        self.setting['margin'] = max(m, max(size)/5.)
      edge = np.array([min(self.molecule.R[:,i])\
        for i in range(3)])
      self.molecule.shift(self.setting['margin']-edge)
      size = 2*self.setting['margin'] + size
      self.setting['celldm'] = np.append(size, [0, 0, 0])
    else:
      self.setting['celldm'] = self.molecule.celldm

