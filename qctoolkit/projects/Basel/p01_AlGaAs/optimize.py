import qctoolkit as qtk
import qctoolkit.ccs as qcs
import qctoolkit.optimization as qop
import copy, shutil, os, glob

def AlGaX_EvOpt(structure, vacancy_ind, ccs_span, **kwargs):
  if 'QMInp' in kwargs:
    inpp = kwargs['QMInp']
  else:
    inpp = qtk.QMInp(structure, 'cpmd')
  if 'T' in kwargs:
    _T = kwargs['T']
  else:
    _T = 300
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

  def clean_file(old_file):
    try:
      os.remove(old_file)
    except:
      pass
  def clean_folder(old_folder):
    try:
      shutil.rmtree(old_folder)
    except:
      pass

  for to_clean in glob.glob('*.inp'):
    clean_file(to_clean)
#  clean_folder('pref')
#  clean_folder('freeAtom')
#  clean_folder('vref')

  ccs = qcs.MoleculeSpan(structure, ccs_span)

  molp = qtk.Molecule()
  inpp.structure = molp.read(structure)
  inpp.setInfo('Ev_per_ref')
  inpp.write('pref.inp')

  inpv = copy.deepcopy(inpp)
  inpv.removeAtom(vacancy_ind)
  inpv.setChargeMultiplicity(inpp.inp.structure.charge, 2)
  inpv.setInfo('Ev_vac_ref')
  inpv.write('vref.inp')

  inpa = copy.deepcopy(inpp)  
  inpa.isolateAtoms(vacancy_ind)
  inpa.setChargeMultiplicity(inpp.inp.structure.charge, 2)
  inpa.setInfo('freeAtom')
  inpa.write('freeAtom.inp')
  
  qtk.QMRun('pref.inp', inpp.program,
            threads=_threads,
            save_restart = True)
  qtk.QMRun('vref.inp', inpp.program,
            threads=_threads,
            save_restart = True)
  if os.path.exists('freeAtom'):
    freeAtomOut = qtk.QMOut('freeAtom/freeAtom.out',inpp.program)
  else:
    freeAtomOut = qtk.QMRun('freeAtom.inp', inpp.program,
                            threads=_threads,
                            save_restart = True,
                            QMReturn=True)

  tmp , init_ccs_coord = ccs.random()
  input_list = [ccs, vacancy_ind, 
                {'QMInp':inpp, 
                 'pref':'pref', 'vref':'vref',
                 'freeAtomE':freeAtomOut.Et, 
                 'threads':_threads}]

  def genCCSInp():
    _tmp, _coord = ccs.random()
    return _coord
 
  mcopt = qop.MonteCarlo(Ev_ccs, input_list, genCCSInp, 
                         power=1, log_file=logfile)
  mcopt.run()

#  qcs.optimize.mc(Ev_ccs, init_ccs_coord, ccs, input_list, 
#                  target=_target, T=_T)
#
#  out = Ev_ccs(init_ccs_coord, ccs, vacancy_ind, QMInp=inpp,
#               pref='pref', vref='vref', 
#               freeAtomE=freeAtomOut.Et,
#               threads = _threads)

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
  base_inp.debug()

  if 'pref' in kwargs and 'vref' in kwargs:
    alchem = True
    perfect_ref = kwargs['pref']
    vacancy_ref = kwargs['vref']
  elif 'pref' not in kwargs and 'vref' not in kwargs:
    alchem = False
  else:
    qtk.exit("either 'pref' or 'vref' is missing in kwargs for"\
             "alchemical runs.")
  if 'freeAtomE' in kwargs:
    freeE = kwargs['freeAtomE']
  else:
    freeE = 0

  if 'threads' in kwargs:
    _threads = kwargs['threads']
  else:
    _threads = 1

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

  qtk.progress("Ev_ccs", "running perfect.inp... ")
  if alchem:
    out_wov = qtk.QMRun('perfect.inp', inp_wov.program,
                        threads=_threads,
                        alchemScan=True,
                        alchemRefPath=perfect_ref,
                        alchemRefPrefix='')
  else:
    out_wov = qtk.QMRun('perfect.inp', inp_wov.program,
                        threads=_threads)
  qtk.done(out_wov.Et)
  qtk.progress("Ev_ccs", "running vacancy.inp... ")
  if alchem:
    out_wv = qtk.QMRun('vacancy.inp', inp_wv.program,
                        threads=_threads,
                        alchemScan=True,
                        alchemRefPath=vacancy_ref,
                        alchemRefPrefix='')
  else:
    out_wv = qtk.QMRun('vacancy.inp', inp_wv.program,
                        threads=_threads)
  qtk.done(out_wv.Et)

  qtk.report("Ev_ccs", "Ev=", out_wov.Et - out_wv.Et - freeE)
  return out_wov.Et - out_wv.Et - freeE
