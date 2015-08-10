import qctoolkit as qtk
import copy, shutil, os

def Ev_ccs(ccs_span, ccs_coord, vacancy_index, **kwargs):
  """
  single point calculation of vacancy energy in crystal
  either reference (true) or predicted (pred) calculations
  can be assigned

  vacancy_index starts from 1
  """
  if 'QMInp' not in kwargs:
    qtk.exit("kwargs: 'QMInp' is missing.\n"\
             + "It should be set to QMInp object of "\
             + "system without vacancies.\n"\
             + "It is necessary for inp settings")
  base_inp = kwargs['QMInp']

  if 'threads' in kwargs:
    _threads = kwargs['threads']
  else:
    _threads = 1

  if 'no_report' in kwargs and kwargs['no_report']:
    _no_report = True
  else:
    _no_report = False

  inp_wov = copy.deepcopy(base_inp)
  inp_wv = copy.deepcopy(base_inp)
  inp_wov.inp.structure = ccs_span.generate(**ccs_coord)
  inp_wv.inp.structure = ccs_span.generate(**ccs_coord)

  inp_wv.removeAtom(vacancy_index)
  inp_wv.setChargeMultiplicity(inp_wov.inp.structure.charge, 2)

  if os.path.exists('perfect'):
    shutil.rmtree('perfect')
  if os.path.exists('vacancy'):
    shutil.rmtree('vacancy')

  inp_wov.write('perfect.inp', no_warning=True)
  inp_wv.write('vacancy.inp', no_warning=True)

  if not _no_report:
    qtk.progress("Ev_ccs", "running perfect.inp")
  out_wov = qtk.QMRun('perfect.inp', inp_wov.program,
                      threads=_threads,
                      QMReturn=True,
                      cleanup=True)
  if not _no_report:
    qtk.done()
    qtk.progress("Ev_ccs", "running vacancy.inp")
  out_wv = qtk.QMRun('vacancy.inp', inp_wv.program,
                      threads=_threads,
                      QMReturn=True,
                      cleanup=True)
  if not _no_report:
    qtk.done()

  return out_wov.Et - out_wv.Et
