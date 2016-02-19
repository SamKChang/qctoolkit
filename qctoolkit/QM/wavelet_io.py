import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput

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

class WaveletOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
