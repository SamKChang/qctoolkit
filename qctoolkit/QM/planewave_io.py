import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
import universal as univ
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

    if 'cutoff' not in kwargs:
      self.setting['cutoff'] = 100
    if not self.setting['periodic'] and 'isolation' not in kwargs:
      self.setting['isolation'] = 'mt'
    self.pp_files = []
    if 'periodic' in self.setting and self.setting['periodic']:
      self.celldm2lattice()
    if 'pp_type' not in kwargs:
      self.setting['pp_type'] = 'geodecker'
    if 'full_kmesh' not in self.setting:
      self.setting['full_kmesh'] = False

    univ.getCelldm(self) 

  def write(self, name=None, **kwargs):
    if self.setting['periodic']:
      self.celldm2lattice()
    inp, molecule = \
      GenericQMInput.write(self, name, **self.setting)

    if 'pp_list' in self.setting:
      pp_list = self.setting['pp_list']
      itr = 1
      for pp_data in pp_list:
        pp_inds = pp_data[0]
        if type(pp_inds) is not list:
          pp_inds = [pp_inds]
        pp_name = pp_data[1]
        pp = pp_data[2]
        pp.setting['program'] = self.setting['program']
        pp.write(pp_name, inplace=False)
        molecule.setAtoms(pp_inds, string=pp_name)
        Zn = molecule.type_list[pp_inds[0]]
        molecule.setAtoms(pp_inds, element=Zn + str(itr))
        itr += 1

    return inp, molecule

  def celldm2lattice(self):
    cd = self.setting['celldm']
    if 'scale' in self.setting:
      sc = self.setting['scale']
    else:
      sc = [1.0 for i in range(3)]
    self.setting['lattice'] = qtk.celldm2lattice(cd, scale=sc)

  def cornerMargin(self, *args, **kwargs):
    pass

class PlanewaveOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
