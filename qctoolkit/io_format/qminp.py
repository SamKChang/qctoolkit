import qctoolkit as qtk
import qctoolkit.io_format.setting_pw as pw

class PwInp(object):
  def __init__(self, structure_inp, info, **kwargs):
    self.setting = pw.Setting()
    self.atom_list = {}
    self.structure = qtk.geometry.Molecule()
    if type(structure_inp) == str:
      self.structure.read(structure_inp, **kwargs)
    else:
      self.structure = copy.deepcopy(structure_inp)
    if self.structure.scale:
      self.setting.scale = self.structure.scale
      self.setting.set_scale = True
      self.setting.isolated = False
    self.setting.celldm = self.structure.celldm
    self.info = info

  def load_structure(self, new_structure, **kwargs):
    self.structure.read(new_structure, **kwargs)
    if self.setting.set_multiplicity:
      _multiplicity = self.setting.multiplicity
    else:
      _multiplicity = 'auto'
    if self.setting.set_charge:
      _charge = self.setting.charge
    else:
      _charge = 'auto'

    print self.setting.set_multiplicity
    print _charge, _multiplicity
    self.structure.setChargeMultiplicity(_charge,
                                         _multiplicity,
                                         **kwargs)

