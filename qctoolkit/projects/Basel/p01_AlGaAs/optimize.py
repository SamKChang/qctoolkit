import qctoolkit as qtk
import copy

def Ev_ccs(ccs_span, ccs_coord, vacancy_index, **kwargs):
  """
  single point calculation of vacancy energy in crystal
  either reference (true) or predicted (pred) calculations
  can be assigned
  """
  print ccs_coord

  if 'QMInp' not in kwargs:
    qtk.exit("kwargs: 'QMInp' is missing.\n"\
             + "It should be set to QMInp object of "\
             + "system without vacancies.\n"\
             + "It is necessary for inp settings")

  inp_wov = copy.deepcopy(kwargs['QMInp'])
  inp_wv = copy.deepcopy(kwargs['QMInp'])
  inp_wov.inp.structure = ccs_span.generate(**ccs_coord)
  inp_wv.inp.structure = inp_wov.inp.structure.remove_atom(vacancy_index)
  print inp_wov.inp.structure.N, inp_wv.inp.structure.N

  inp_wov.write('stdout')
  inp_wv.write('stdout')

  return inp_wov, inp_wv
