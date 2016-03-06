import qctoolkit as qtk
import qctoolkit.ccs as qcs
import qctoolkit.optimization as qop
import qctoolkit.alchemy as qal
import copy, shutil, os, glob

def AlGaX_EvOpt(structure, vacancy_ind, ccs_span, **kwargs):

  qm_setting = {}
  if 'qm_setting' in kwargs:
    qm_setting = kwargs['qm_setting']
  qm_setting['save_restart'] = True

  if 'QMInp' in kwargs:
    baseinp = kwargs['QMInp']
  else:
    baseinp = qtk.QMInp(structure, program='cpmd')
  if 'T' in kwargs:
    _T = kwargs['T']
  else:
    _T = 1
  if 'target' in kwargs:
    _target = kwargs['target']
  else:
    _target = 0

  if 'log_file' in kwargs:
    logfile = kwargs['log_file']
  else:
    logfile = 'AlGaX_EvOpt.log'

  if 'threads' in kwargs:
    _threads = kwargs['threads']
  else:
    _threads = qtk.cpu_count

  if 'threads_per_job' in kwargs:
    _threadspj = kwargs['threads_per_job']
  else:
    _threadspj = _threads
  _parallel = int(_threads/_threadspj)

  if 'optimizer' in kwargs:
    _optimizer = kwargs['optimizer']
    if _optimizer == 'GA':
      if 'population_size' in kwargs:
        _population_size = kwargs['population_size']
      else:
        _population_size = qtk.setting.cpu_count
  else:
    _optimizer = 'MC'

  ccs = qtk.CCS(structure, ccs_span)

  inpp = qtk.QMInp(structure, **qm_setting)
  inpp.setting['info'] = 'Ev_per_ref'
  if not os.path.exists('pref/pref.out'):
    inpp.run('pref')

  inpv = qtk.QMInp(structure, **qm_setting)
  inpv.removeAtoms(vacancy_ind)
  inpv.setChargeMultiplicity(0, 2)
  inpv.setting['info'] = 'Ev_vac_ref'
  if not os.path.exists('vref/vref.out'):
    inpv.run('vref')

  inpa = qtk.QMInp(structure, **qm_setting)
  inpa.isolateAtoms(vacancy_ind)
  inpa.setChargeMultiplicity(0, 2)
  inpa.setting['info'] = 'freeAtom'
  if not os.path.exists('freeAtom/freeAtom.out'):
    inpa.run('freeAtom')
  freeAtomOut = qtk.QMOut('freeAtom/freeAtom.out')

  tmp , init_ccs_coord = ccs.random()

  qm_setting['threads'] = _threadspj
  penalty_setting = {
                      'QMInp':baseinp,
                      'freeAtomE':freeAtomOut.Et, 
                      'qm_setting': qm_setting,
                    }
  if 'alchemy' in kwargs and kwargs['alchemy']:
    penalty_setting['pref'] = 'pref'
    penalty_setting['vref'] = 'vref'
  input_list = [ccs, vacancy_ind, penalty_setting]
  
  def genCCSInp():
    _coord = ccs.random()[1]
    return _coord
 
  op_setting = {
                 'power': 1, 
                 'log_file': logfile,
                 'target': _target,
                 'parallel': _parallel,
                 'T': _T,
               }

  qtk.report('start optimizer')
  if _optimizer == 'MC':
    cylopt = qop.MonteCarlo(Ev_ccs, input_list, genCCSInp, 
                            **op_setting)
  elif _optimizer == 'GA':
    cylopt = qop.GeneticOptimizer(Ev_ccs, input_list, genCCSInp, 
                                  ccs.mate, _population_size,
                                  **op_setting)
  qtk.report('optimizer initialized')

  cylopt.run()

def Ev_ccs(ccs_coord, ccs_span, vacancy_index, **kwargs):
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

  qm_setting = {}
  if 'qm_setting' in kwargs:
    qm_setting = kwargs['qm_setting']

  if 'pref' in kwargs and 'vref' in kwargs:
    alchem = True
    perfect_ref = kwargs['pref']
    vacancy_ref = kwargs['vref']
  elif 'pref' not in kwargs and 'vref' not in kwargs:
    alchem = False

  freeE = qtk.QMOut('freeAtom/freeAtom.out')
  freeE.inUnit('ev')

  if 'threads' in kwargs:
    _threads = kwargs['threads']
  else:
    _threads = 1

  inp_wov = qtk.QMInp(ccs_span.generate(**ccs_coord))
  inp_wv = qtk.QMInp(ccs_span.generate(**ccs_coord))

  inp_wv.removeAtoms(vacancy_index)
  inp_wv.setChargeMultiplicity(0, 2)

  perfect = 'ev_perfect' + str(os.getpid())
  vacancy = 'ev_vacancy' + str(os.getpid())
  perfectinp = perfect + '.inp'
  vacancyinp = vacancy + '.inp'
  inp_wov.molecule.name = perfectinp
  inp_wv.molecule.name = vacancyinp

  if os.path.exists(perfect):
    shutil.rmtree(perfect)
  if os.path.exists(vacancy):
    shutil.rmtree(vacancy)

  print ccs_coord
  if alchem:
    out_wov = qtk.Al1st(inp_wov, ref_dir=perfect_ref, **qm_setting)
    out_wv = qtk.Al1st(inp_wv, ref_dir=vacancy_ref, **qm_setting)
  else:
    out_wov = inp_wov.run(**qm_setting)
    out_wv = inp_wv.run(**qm_setting)
  try:
    os.remove(perfectinp)
    os.remove(vacancyinp)
  except OSError:
    shutil.rmtree(perfectinp)
    shutil.rmtree(vacancyinp)

  out_wov.inUnit('ev')
  out_wv.inUnit('ev')

  final = out_wov - out_wv - freeE

  msg = str(out_wov.Et) + '-(' + str(out_wv.Et) + \
        '+' + str(freeE.Et) + ') = ' + str(final.Et)
  qtk.report('trial Ev', msg)

  return final.Et
