import qctoolkit as qtk
import copy
from general_io import GenericMDInput
from general_io import GenericMDOutput

class AIMDInp(GenericMDInput):
  def __init__(self, molecule, **kwargs):
    GenericMDInput.__init__(molecule, **kwargs)
    if 'theory' not in kwargs: 
      self.setting['theory'] = 'PBE'

    if 'vdw' not in kwargs: 
      self.setting['vdw'] = kwargs['vdw']

    if 'cutoff' not in kwargs:
      self.setting['cutoff'] = 100

    if 'periodic' not in kwargs: 
      self.setting['periodic'] = True

    if 'em_step' not in kwargs: 
      self.setting['em_step'] = 200

    if 'eq_step' not in kwargs: 
      self.setting['eq_step'] = 200

    if 'md_step' not in kwargs: 
      self.setting['md_step'] = 5000

    if 'md_mode' not in kwargs: 
      self.setting['md_mode'] = 'BOMD'

    # celldm can always be overwritten
    if 'celldm' not in kwargs:
      # for the case of molecule xyz input (molecule)
      # set default orthorombic celldm
      if not self.molecule.celldm:
        box = self.molecule.getBox()
        # set defualt margin in specified in setting.py, 
        # grow with box size
        if 'margin' not in kwargs:
          m = qtk.setting.box_margin
          self.setting['margin'] = max(m, max(box)/5.)
        edge = np.array([min(self.molecule.R[:,i])\
          for i in range(3)])
        self.molecule.shift(self.setting['margin']-edge)
        box = 2*self.setting['margin'] + box
        self.setting['celldm'] = np.append(box, [0, 0, 0])
      else:
        self.setting['celldm'] = self.molecule.celldm


    for key in self.md_setting.iterkeys():
      qtk.report("MDJob", "setting", key, 
                 ":", self.md_setting[key])

    em_setting = copy.deepcopy(self.setting)
    em_setting['md_step'] = em_setting['em_step']
    em_setting['mode'] = 'geopt'
    em_setting['info'] = ' MD em job with ' + molecule
    self.emInp = qtk.QMInp(molecule, **em_setting)

    eq_setting = copy.deepcopy(self.setting)
    eq_setting['md_step'] = eq_setting['eq_step']
    eq_setting['info'] = ' MD eq job with ' + molecule
    eq_setting['mode'] = self.md_setting['md_mode']
    self.eqInp = qtk.QMInp(molecule, **eq_setting)

    md_setting = copy.deepcopy(self.setting)
    md_setting['mode'] = self.md_setting['md_mode']
    md_setting['info'] = ' MD md job with ' + molecule
    self.mdInp = qtk.QMInp(molecule, **md_setting)

class AIMDOut(GenericMDOutput):
  def __init__(self, out_dir, **kwargs):
    GenericMDOutput.__init__(self, out_dir, **kwargs)
