import qctoolkit as qtk
import numpy as np

def getCelldm(self):
  if 'celldm' not in self.setting:
    if not self.molecule.celldm:
      box = self.molecule.getBox()
      if 'margin' not in self.setting:
        m = qtk.setting.pw_margin
        self.setting['margin'] = max(m, max(box)/5.)
      edge = np.array([min(self.molecule.R[:,i])\
        for i in range(3)])
      self.molecule.shift(self.setting['margin']-edge)
      box = 2*self.setting['margin'] + box
      self.setting['celldm'] = np.append(box, [0, 0, 0])
    else:
      self.setting['celldm'] = self.molecule.celldm
