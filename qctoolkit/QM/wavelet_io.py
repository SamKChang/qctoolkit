import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
import universal as univ
import copy

class WaveletInput(GenericQMInput):
  """
  From WaveletBasis Input:
  generic class holder for wavelet qmcode. It provide basic
  default settings.
  ===
  """
  __doc__ = GenericQMInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if self.setting['periodic']:
      univ.getCelldm(self)
      box = copy.deepcopy(self.setting['celldm'][:3])
      self.setting['box'] = box
    else:
      self.setting['box'] = False
      self.setting['celldm'] = False

class WaveletOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
