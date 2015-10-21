import qctoolkit as qtk
import copy

class MDInp(object):
  def __init__(self, molecule, **kwargs):
    self.md_setting = {}
    if 'theory' in kwargs: 
      self.md_setting['theory'] = kwargs['theory']
    else: 
      self.md_setting['theory'] = 'PBE'

    if 'vdw' in kwargs: 
      self.md_setting['vdw'] = kwargs['vdw']
    else: 
      self.md_setting['vdw'] = False

    if 'cutoff' in kwargs:
      self.md_setting['cutoff'] = kwargs['cutoff']
    else:
      self.md_setting['cutoff'] = 100

    if 'celldm' in kwargs: 
      self.md_setting['celldm'] = kwargs['celldm']
    else: 
      self.md_setting['celldm'] = [20, 20, 20, 0, 0, 0]

    if 'periodic' in kwargs: 
      self.md_setting['periodic'] = kwargs['periodic']
    else: 
      self.md_setting['periodic'] = True

    if 'em_step' in kwargs: 
      self.md_setting['em_step'] = kwargs['em_step']
    else: 
      self.md_setting['em_step'] = 200

    if 'eq_step' in kwargs: 
      self.md_setting['eq_step'] = kwargs['eq_step']
    else: 
      self.md_setting['eq_step'] = 200

    if 'md_step' in kwargs: 
      self.md_setting['md_step'] = kwargs['md_step']
    else: 
      self.md_setting['md_step'] = 5000

    if 'temperature' in kwargs: 
      self.md_setting['temperature'] = kwargs['temperature']
    else: 
      self.md_setting['temperature'] = 298

    if 'temperature_tolerance' in kwargs:
      self.md_setting['temperature_tolerance'] = \
      kwargs['temperature_tolerance']
    else: 
      self.md_setting['temperature_tolerance'] = 100

    if 'md_mode' in kwargs: 
      self.md_setting['md_mode'] = kwargs['md_mode']
    else: 
      self.md_setting['md_mode'] = 'BOMD'

    if 'md_sample_rate' in kwargs:
      self.md_setting['md_sample_rate'] = \
      kwargs['md_sample_rate']
    if 'theory' in kwargs:
      self.md_setting['md_sample_rate'] = 10

    for key in self.md_setting.iterkeys():
      qtk.report("MDJob", "setting", key, 
                 ":", self.md_setting[key])

    em_setting = copy.deepcopy(self.md_setting)
    em_setting['md_step'] = em_setting['em_step']
    em_setting['mode'] = 'geopt'
    em_setting['info'] = ' MD em job with ' + molecule
    self.emInp = qtk.QMInp(molecule, **em_setting)

    eq_setting = copy.deepcopy(self.md_setting)
    eq_setting['md_step'] = eq_setting['eq_step']
    eq_setting['info'] = ' MD eq job with ' + molecule
    eq_setting['mode'] = self.md_setting['md_mode']
    self.eqInp = qtk.QMInp(molecule, **eq_setting)

    md_setting = copy.deepcopy(self.md_setting)
    md_setting['mode'] = self.md_setting['md_mode']
    md_setting['info'] = ' MD md job with ' + molecule
    self.mdInp = qtk.QMInp(molecule, **md_setting)
