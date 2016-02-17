import qctoolkit as qtk
import qctoolkit.QM.qmcode.cpmd as cpmd
from qctoolkit.QM.qmdir import qmDir
import subprocess as sp
import os, re, shutil
import alpath as alp

def ReferenceRun(inp, program=qtk.setting.qmcode, **kwargs):
  if program == 'cpmd':
    # to pass 'save_restart' explicitly
    if 'save_restart' in kwargs:
      del kwargs['save_restart']
    qtk.QMRun(inp, 'cpmd', 
              save_restart=True, 
              **kwargs)

def FirstOrderRun(inp, program=qtk.setting.qmcode, **kwargs):
  if program == 'cpmd':
    inpdir, inpname, psinp, new_run, kwargs\
      = qmDir(inp, **kwargs)
    if new_run:   
      if 'rst' in kwargs:
        rst = kwargs['rst']
        del kwargs['rst']
      else:
        rst = 'RESTART.1'
      ref = kwargs['ref_path']
      rst_src = os.path.join(ref, rst)
      if not os.path.exists(ref):
        qtk.exit("FirstOrderRun: ref_path", ref, "not found")
      if not os.path.exists(rst_src):
        qtk.exit("FirstOrderRun: RESTART file", 
                 rst_src, "not found")
      rst_trg = os.path.join(inpdir, 'RESTART')
      os.link(rst_src, rst_trg)
  
      cwd = os.getcwd()
      qout = qtk.QMRun(inpname, 'cpmd', inplace=True, 
                       restart=True, maxstep=1, **kwargs)
      return qout

# not working
def SecondOrderRun(inp, program=qtk.setting.qmcode, **kwargs):
  if program == 'cpmd':
    inpdir, inpname, psinp, new_run, kwargs\
      = qmDir(inp, **kwargs)
    if new_run:
      if 'inp_base' in kwargs:
        _inp_base = kwargs['inp_base']
      else:
        _inp_base = None
      if 'dl' in kwargs:
        dl = kwargs['dl']
      else:
        dl = 0.05
      dln = '_%03d' % (dl*100)
      inpname_dl = inpdir + '/' + psinp +   dln  + '.inp'
      inpname_L  = inpdir + '/' + psinp + '_100' + '.inp'
      ref_i  = inpdir + '/' + psinp + '_000-d0.inp'
      ref_dl = inpdir + '/' + psinp + '_000-dl.inp'
      tar_i  = inpdir + '/' + psinp + '_100-d0.inp'
      tar_dl = inpdir + '/' + psinp + '_100-dl.inp'
      if 'rst' in kwargs:
        rst = kwargs['rst']
        del kwargs['rst']
      else:
        rst = 'RESTART.1'
      ref = kwargs['ref_path']
      if 'ref_inp' in kwargs:
        rip = kwargs['ref_inp']
      else:
        refroot = re.sub('.*/','',ref)
        rip = ref + '/' + refroot + '.inp'
      rst_src = ref + '/' + rst
      if not os.path.exists(ref):
        qtk.exit("SecondOrderRun: ref_path="+ref+" not found")
      if not os.path.exists(rst_src):
        qtk.exit("SecondOrderRun: RESTART file="+\
                 rst_src+" not found")
      if not os.path.exists(rip):
        qtk.exit("SecondOrderRun: ref_inp="+rip+" not found")
      rst_trg = inpdir + '/RESTART_d0'
      os.link(rst_src, rst_trg)

      # prepare input file of dl full SCF
      if qtk.matching(inpname, r'.*_0*\.psp'):
        shutil.copy(inpname, inpname_dl)
        shutil.copy(inpname, inpname_L)
        qtk.replace(inpname_dl, '_0*', dln)
        qtk.replace(inpname_L , '_0*', '_100')
      else:
        xyz_i = qtk.Molecule()
        xyz_f = qtk.Molecule()
        xyz_i.read_cpmdinp(rip)
        xyz_f.read_cpmdinp(inpname)
        alpath = alp.PathScan(xyz_i, xyz_f, 'cpmd', 
                              inp_base=_inp_base)
        inp_dl1 = alpath.l(dl)
        inp_L   = alpath.l(1)
        inp_dl1.setCutoff(20)
        inp_L.setCutoff(20)
        inp_dl1.write(inpname_dl)
        inp_L.write(inpname_L)
      os.link(inpdir+'/RESTART_d0', inpdir+'/RESTART')
      qdl = qtk.QMRun(inpname_dl, 'cpmd', inplace=True, 
                      restart=True, save_restart=True, **kwargs)
      os.rename(inpdir+'/RESTART.1', inpdir+'/RESTART_dl')

      # restart at lambda=0
      os.rename(inpname_L, tar_i)
      qti  = qtk.QMRun(tar_i,  'cpmd', inplace=True,
                       save_restart=True,
                       maxstep=1, restart=True, **kwargs)
#      os.rename(inpname  , ref_i)
#      qri  = qtk.QMRun(ref_i,  'cpmd', inplace=True, 
#                       save_restart=True,
#                       maxstep=1, restart=True, **kwargs)

      # restart at lambda=dl
      os.remove(inpdir+'/RESTART')
      shutil.copy(tar_i, tar_dl)
      os.link(inpdir+'/RESTART_dl', inpdir+'/RESTART')
      qtd  = qtk.QMRun(tar_dl,  'cpmd', inplace=True, 
                       save_restart=True,
                       maxstep=1, restart=True, **kwargs)
      shutil.copy(ref_i, ref_dl)
      qrd  = qtk.QMRun(ref_dl,  'cpmd', inplace=True, 
                       save_restart=True,
                       maxstep=1, restart=True, **kwargs)
      # clean up
      os.remove(inpdir+'/RESTART')
      os.remove(inpdir+'/RESTART_dl')
      os.remove(inpdir+'/RESTART')
