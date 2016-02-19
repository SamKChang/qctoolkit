import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
import numpy as np

class PlanewaveInput(GenericQMInput):
  """
  From PlanwaveInput:
  generic class holder for plane wave qmcode. It provide basic
  default settings.
  """
  __doc__ = GenericQMInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if 'periodic' not in kwargs:
      self.setting['periodic'] = True
    if 'symmetry' not in kwargs and self.setting['periodic']:
      self.setting['symmetry'] = 'orthorhombic'
    if 'cutoff' not in kwargs:
      self.setting['cutoff'] = 100
    if not self.setting['periodic'] and 'isolation' not in kwargs:
      self.setting['isolation'] = 'mt'
    self.pp_files = []
    self.pp_path = ''
  
    # celldm can always be overwritten
    if 'celldm' not in kwargs:
      # for the case of molecule xyz input (molecule)
      # set default orthorombic celldm
      if not self.molecule.celldm:
        box = self.molecule.getBox()
        # set defualt margin in specified in setting.py, 
        # grow with box size
        if 'margin' not in kwargs:
          m = qtk.setting.pw_margin
          self.setting['margin'] = max(m, max(box)/5.)
        edge = np.array([min(self.molecule.R[:,i])\
          for i in range(3)])
        self.molecule.shift(self.setting['margin']-edge)
        box = 2*self.setting['margin'] + box
        self.setting['celldm'] = np.append(box, [0, 0, 0])
      else:
        self.setting['celldm'] = self.molecule.celldm

  def celldm2lattice(self):
    cd = self.setting['celldm']
    if cd[3]==0 and cd[4]==0 and cd[5]==0:
      self.setting['lattice'] = np.array([[cd[0],   0.0,   0.0],
                                          [  0.0, cd[1],   0.0],
                                          [  0.0,   0.0, cd[2]]])
  def center(self, coord):
    self.molecule.center(coord)
  def shift(self, coord):
    self.molecule.shift(coord)

class PlanewaveOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
