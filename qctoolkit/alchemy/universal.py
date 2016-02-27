import qctoolkit as qtk

def toInp(qminp, **kwargs):
  inpClass = qtk.QM.general_io.GenericQMInput
  is_qminp = issubclass(qminp.__class__, inpClass)
  if not is_qminp:
    if type(qminp) is str or type(qminp) is qtk.Molecule:
      qminp = qtk.QMInp(qminp, **kwargs)
  return qminp
