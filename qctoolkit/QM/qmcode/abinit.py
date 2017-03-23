import qctoolkit as qtk
from qctoolkit.QM.planewave_io import PlanewaveInput
from qctoolkit.QM.planewave_io import PlanewaveOutput
import os, copy, glob, re
import qctoolkit.QM.qmjob as qmjob
import numpy as np
import universal as univ
from bigdft import PPCheck
from collections import OrderedDict as odict
from numbers import Number
import subprocess as sp

class inp(PlanewaveInput):
  """
  abinit input class.
  """
  __doc__ = PlanewaveInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    PlanewaveInput.__init__(self, molecule, **kwargs)
    self.setting.update(**kwargs)
    self.backup()
    self.content = odict()
    self.content['datasets'] = odict([('ndtset', 1)])

  def write(self, name=None, **kwargs):

    def writeInp(name=None, **setting):
      self.setting.update(setting)
      self.setting['no_molecule'] = False
      inp, molecule = \
        super(PlanewaveInput, self).write(name, **self.setting)

      molecule.sort()
      type_index = molecule.index
      type_list = molecule.type_list

      self.pp_path = qtk.setting.bigdft_pp
      if 'pp_path' in self.setting:
        self.pp_path = self.setting['pp_path']

      if 'pp_theory' not in self.setting:
        self.setting['pp_theory'] = self.setting['theory']

      pp_files = []
      typat = []
      for a in range(len(type_index)-1):
        start = type_index[a]
        end = type_index[a+1]
        for i in range(end - start):
          typat.append(a+1)
      for i in range(molecule.N):
        pp_file = 'psppar.' + molecule.type_list[i]
        pp_list = set([pp[1] for pp in pp_files])
        if pp_file not in pp_list:
          pp_src = PPCheck(self.setting['theory'],
                           self.setting['pp_theory'],
                           self.pp_path,
                           molecule.type_list[i])
          pp_files.append([pp_src, pp_file])

      ##################
      # scf section #
      ##################
      self.content['scf'] = odict()

      # restart check #
      if 'restart' in self.setting and self.setting['restart']\
      and 'band_scan' not in self.setting\
      and 'dos_mesh' not in self.setting:
        self.content['scf']['irdwfk'] = 1
        self.content['scf']['getwfk'] = -1
      if 'restart_density' in self.setting \
      and self.setting['restart_density']\
      and 'band_scan' not in self.setting\
      and 'dos_mesh' not in self.setting:
        self.content['scf']['irdden'] = 1
        self.content['scf']['getden'] = -1

      if 'kmesh' not in self.setting:
        self.setting['kmesh'] = [1,1,1]
      if 'kshift' not in self.setting:
        self.setting['kshift'] = [0.0,0.0,0.0]

      # kpoints check #
      if 'kpoints' in self.setting:
        self.content['scf']['kptopt'] = 0
        self.content['scf']['nkpt'] = len(self.setting['kpoints'])
        self.content['scf']['kpt'] = self.setting['kpoints']
      else:
        if self.setting['full_kmesh']:
          self.content['scf']['kptopt'] = 3
        self.content['scf']['ngkpt'] = self.setting['kmesh']
        if len(np.array(self.setting['kshift']).shape) > 1:
          self.content['scf']['nshiftk'] = len(self.setting['kshift'])
        self.content['scf']['shiftk'] = self.setting['kshift']
        nbnd = None
      if 'ks_states' in self.setting and self.setting['ks_states']:
        vs = int(round(molecule.getValenceElectrons() / 2.0))
        nbnd = self.setting['ks_states'] + vs
        if 'd_shell' in self.setting:
          for a in molecule.type_list:
            if a in self.setting['d_shell'] and qtk.n2ve(a) < 10:
              nbnd += 5
        if 'band_scan' not in self.setting \
        and 'dos_mesh' not in self.setting:
          self.content['scf']['nband'] = nbnd

      # system setup #
      if molecule.charge != 0:
        self.content['scf']['charge='] = molecule.charge
      if molecule.multiplicity != 1:
        self.content['scf']['nsppol'] = 2
        self.content['scf']['occopt'] = 7
      if 'save_restart' not in self.setting \
      or not self.setting['save_restart']:
        if 'restart' in self.setting and self.setting['restart']\
        and ('band_scan' not in self.setting\
        and 'dos_mesh' not in self.setting):
          self.content['scf']['prtwf'] = 0
      if 'wf_convergence' in self.setting:
        if 'restart' in self.setting and self.setting['restart']\
        and ('band_scan' not in self.setting\
        and 'dos_mesh' not in self.setting):
          self.content['scf']['tolwfr'] = \
          self.setting['wf_convergence']

      # clean up for the case of restart
      if 'restart' in self.setting and self.setting['restart']\
      and ('band_scan' in self.setting or 'dos_mesh' in self.setting):
        keep_lst = [
          'nband',
          'tolwfr',
          'nstep',
          'ecut',
          'irdwfk',
          'irdden',
          'getwfk',
          'getden',
        ]
        for key in self.content['scf'].iterkeys():
          if key not in keep_lst and 'prt' not in key:
            self.content['scf'].pop(key)

      ######################
      # additional setting #
      ######################
      if 'abinit_setting' in self.setting:
        self.content['additional'] = odict([])
        for item in self.setting['abinit_setting']:
          self.content['additional'][item] = ' '

      def is_strnum(q):
        if isinstance(q, Number) or type(q) is str:
          return True
        else:
          return False

      #################
      # bandstructure #
      #################
      if 'band_scan' in self.setting and self.setting['band_scan']:

        if molecule.symmetry and molecule.symmetry.lower() == 'fcc':
          if 'kshift' not in self.setting:
            self.setting['kshift'] = [
              [0.5, 0.5, 0.5],
              [0.5, 0.0, 0.0],
              [0.0, 0.5, 0.0],
              [0.0, 0.0, 0.5],
            ]

        bnds = self.setting['band_scan']
        if len(bnds) != 2 \
        or is_strnum(bnds[0]) \
        or len(bnds[0]) != len(bnds[1]) - 1:
          qtk.exit('band_scan format: [lst_div, lst_coord]')

        bnd_content = odict([
          ('iscf', -2),
          ('getden', -1),
          ('kptopt', -len(bnds[0])),
          ('tolwfr', self.setting['wf_convergence']),
          ('ndivk', bnds[0]),
          ('kptbounds', bnds[1]),
        ])
        if nbnd:
          bnd_content['nband'] = nbnd
        if 'save_restart' not in self.setting \
        or not self.setting['save_restart']:
          bnd_content['prtwf'] = 0
        else:
          bnd_content['prtwf'] = 1
        if 'save_density' not in self.setting \
        or not self.setting['save_density']:
          bnd_content['prtden'] = 0
        else:
          bnd_content['prtden'] = 1
        if 'dos_mesh' in self.setting:
          bnd_content['prtden'] = 1

        if 'restart' in self.setting and self.setting['restart']:
          bnd_content['irdden'] = 1
          bnd_content['irdwfk'] = 1
          bnd_content['getwfk'] = -1


        self.content['band_scan'] = bnd_content

      #####################
      # density of states #
      #####################
      if 'dos_mesh' in self.setting and self.setting['dos_mesh']:
        dos_content = odict([
          ('iscf', -3),
          ('ngkpt', self.setting['dos_mesh']),
          ('shiftk', [0.0, 0.0, 0.0]),
          ('tolwfr', self.setting['wf_convergence']),
          ('prtwf', 0),
        ])
        if 'band_scan' in self.setting:
          dos_content['getden'] = -2
        else:
          dos_content['getden'] = -1

        if 'smearing' in self.setting:
          smr = self.setting['smearing']
          smr_dict = {
            'fermi': 3,
            'marzari': 4,
            'marzari_monotonic': 5,
            'paxton': 6,
            'gaussian': 7,
          }
          if smr in smr_dict:
            dos_content['occopt'] = smr_dict[smr]
            dos_content['prtdos'] = 1
          else:
            qtk.warning("%s for occupation scheme not found" % smr)

        if nbnd:
          dos_content['nband'] = nbnd

        if 'save_restart' not in self.setting \
        or not self.setting['save_restart']:
          dos_content['prtwf'] = 0
        else:
          dos_content['prtwf'] = 1
        if 'save_density' not in self.setting \
        or not self.setting['save_density']:
          dos_content['prtden'] = 0
        else:
          dos_content['prtden'] = 1

        if 'restart' in self.setting and self.setting['restart']:
          dos_content['irdden'] = 1
          dos_content['irdwfk'] = 1
          dos_content['getwfk'] = -1

        self.content['dos'] = dos_content

        
      ################
      # cell section #
      ################
      self.content['cell'] = odict([
        ('ecut', float(self.setting['cutoff']/2.)),
        ('nstep', self.setting['scf_step']),
      ])

      if not molecule.symmetry:
        self.content['cell']['acell'] = '3*1.889726124993'
        if 'lattice' not in self.setting:
          self.celldm2lattice()
        lattice_vec = self.setting['lattice']
        self.content['cell']['chkprim'] = 0
      else:
        self.content['cell']['acell'] = [1, 1, 1]
        a0 = molecule.celldm[0] * 1.889726124993
        if  molecule.symmetry.lower() == 'fcc':
          lattice_vec = 0.5 * a0 * np.array([
            [0, 1, 1],
            [1, 0, 1],
            [1, 1, 0],
          ])
        elif  molecule.symmetry.lower() == 'bcc':
          lattice_vec = 0.5 * a0 * np.array([
            [-1, 1, 1],
            [ 1,-1, 1],
            [ 1, 1,-1],
          ])
      self.content['cell']['rprim'] = lattice_vec


      #################
      # atoms section #
      #################
      znucl = []
      for a in range(len(type_index)-1):
        symb = type_list[type_index[a]]
        znucl.append(int(qtk.n2Z(symb)))

      self.content['atoms'] = odict([
        ('natom', molecule.N),
        ('ntypat', (len(type_index) - 1)),
        ('typat', typat),
        ('znucl', znucl),
      ])
      if hasattr(molecule, 'scale') and molecule.scale:
        if hasattr(molecule, 'R_scale'):
          R_scale = molecule.R_scale.copy()
          for i in range(3):
            s = molecule.scale[i]
            R_scale[:,i] = R_scale[:,i] / s
          self.content['atoms']['xred'] = R_scale
        else:
          qtk.warning('R_scale not found but scale is set')
          self.content['atoms']['xangst'] = molecule.R
      else:
        self.content['atoms']['xangst'] = molecule.R
     

      #########################
      # write content to file #
      #########################

      datasets = 1
      for key in ['band_scan', 'dos']:
        if key in self.content:
          datasets += 1

      if self.setting['restart'] and datasets > 1:
        datasets -= 1

      self.content['datasets']['ndtset'] = datasets

      inp.write('# abinit input generated by qctoolkit\n')
      ds_itr = 0
      ds_str = ''

      if self.content['datasets']['ndtset'] == 1:
        self.content['datasets'].pop('ndtset')

      for section_key in self.content.iterkeys():
        if len(self.content[section_key]) > 0:
          inp.write('\n# %s section\n' % section_key)
        if section_key in ['scf', 'band_scan', 'dos'] and datasets > 1:
          ds_itr += 1
          ds_str = str(ds_itr)
        for key, value in self.content[section_key].iteritems():
          if ds_str and section_key in ['scf', 'band_scan', 'dos']:
            key = key + ds_str
          if is_strnum(value):
            inp.write('%s %s\n' % (key, str(value)))
          else:
            inp.write('%s' % key)
            head = ''.join([' ' for _ in key])
            if 'xangst' in key\
            or 'xcart' in key\
            or 'xred' in key:
              inp.write('\n\n')
              for vec in value:
                inp.write('%s' % head)
                for v in vec:
                  inp.write(' %s' % str(v))
                inp.write('\n')
            elif is_strnum(value[0]):
              for v in value:
                inp.write(' %s' % str(v))
              inp.write('\n')
            else:
              for v in value[0]:
                inp.write(' %s' % str(v))
              inp.write('\n')
              for vec in value[1:]:
                inp.write('%s' % head)
                for v in vec:
                  inp.write(' %s' % str(v))
                inp.write('\n')

      if 'restart_wavefunction_path' in self.setting:
        rstwf_path = self.setting['restart_wavefunction_path']
        if os.path.exists(rstwf_path):
          wf_lst = sorted(
            glob.glob(os.path.join(rstwf_path, '*o_*WFK'))
          )
          assert len(wf_lst) > 0
          wf_src = os.path.abspath(wf_lst[-1])
          pp_files.append([wf_src, name + 'i_WFK'])
        else:
          qtk.warning('%s not found' % rstwf_path)

      if 'restart_wavefunction_file' in self.setting:
        rst_file = self.setting['restart_wavefunction_file']
        if os.path.exists(rst_file):
          pp_files.append([rst_file, name + 'i_WFK'])
        else:
          qtk.warning('%s not found' % rst_file)
          
      if 'restart_density_path' in self.setting:
        rstdn_path = self.setting['restart_density_path']
        if os.path.exists(rstdn_path):
          dn_lst = sorted(
            glob.glob(os.path.join(rstdn_path, '*o_*DEN'))
          )
          assert len(dn_lst) > 0
          dn_src = os.path.abspath(dn_lst[-1])
          pp_files.append([dn_src, name + 'i_DEN'])
        else:
          qtk.warning('%s not found' % rstdn_path)

      if 'restart_density_file' in self.setting:
        rst_file = self.setting['restart_density_file']
        if os.path.exists(rst_file):
          pp_files.append([rst_file, name + 'i_DEN'])
        else:
          qtk.warning('%s not found' % rst_file)
          

      inp.close(dependent_files=pp_files)

      if hasattr(inp, 'final_name'):
        self.setting['no_molecule'] = True
        self.setting['root_dir'] = name
        files = \
          super(PlanewaveInput, self).write(name, **self.setting)
        files.extension = 'files'
        files.write(inp.final_name + '\n')
        root = os.path.splitext(inp.final_name)[0]
        files.write(root + '.out\n')
        files.write(root + 'i\n')
        files.write(root + 'o\n')
        files.write(root + 'x\n')
        for pp in pp_files:
          files.write(pp[1] + '\n')
        files.close(no_cleanup = True)

        

      return inp

    setting = copy.deepcopy(self.setting)
    setting.update(kwargs)
    inp = writeInp(name, **setting)
    return inp

  def run(self, name=None, **kwargs):
    self.setting.update(kwargs)
    return univ.runCode(self, PlanewaveInput, name, **self.setting)

class out(PlanewaveOutput):
  def __init__(self, qmout, **kwargs):
    PlanewaveOutput.__init__(self, qmout, **kwargs)
    out_file = open(qmout)
    data = out_file.readlines()
    out_file.close()

    EStrList = filter(lambda x: 'Etotal' in x, data)
    EList = filter(lambda x: 'ETOT' in x, data)
    self.scf_step = len(EList)
    if self.scf_step > 0:
      Et_list = [float(filter(None, s.split(' '))[2]) for s in EList]
      self.Et = Et_list[-1]
    elif len(EStrList) > 0:
      EStr = EStrList[-1]
      self.Et = float(EStr.split(' ')[-1])

    if len(EStrList) > 0:
      EStr = EStrList[-1]
      detailInd = data.index(EStr)
      self.detail = data[detailInd-7:detailInd]

    GStr = filter(lambda x: 'G(1)' in x, data)[-1]
    ind_g = len(data) - data[::-1].index(GStr) - 1
    G_lattice = []
    for i in range(3):
      G_lattice.append(
        np.fromstring(data[ind_g + i].split('=')[-1], sep=' ')
      )
    self.G_lattice = np.array(G_lattice)

    xangst = filter(lambda x: 'xangst' in x, data)[-1]
    angst_n = len(data) - data[::-1].index(xangst) - 1
    xcart = filter(lambda x: 'xcart' in x, data)[-1]
    cart_n = len(data) - data[::-1].index(xcart) - 1
    Rstr = copy.deepcopy(data[angst_n:cart_n])
    Rstr[0] = Rstr[0].replace('xangst', '')
    R = [[float(r) for r in filter(None, s.split(' '))] for s in Rstr]
    N = len(R)
    ZstrOriginal = filter(lambda x: ' typat' in x, data)[-1]
    Zstr = ZstrOriginal.replace('typat', '')
    Zind = [int(z) for z in filter(None, Zstr.split(' '))]
    ZindItr = data.index(ZstrOriginal)
    while len(Zind) != N:
      ZindItr += 1
      ZindNewStr = filter(None, data[ZindItr].split(' '))
      ZindNew = [int(z) for z in ZindNewStr]
      Zind.extend(ZindNew)
    NZnuc = filter(lambda x: 'ntypat' in x, data)[-1]
    NZnuc = int(filter(None, NZnuc.split(' '))[-1])
    Znuc = filter(lambda x: 'znucl ' in x, data)[-1]
    line_znuc = len(data) - data[::-1].index(Znuc) 
    Znuc = filter(None, Znuc.replace('znucl', '').split(' '))
    Znuc = [float(z) for z in Znuc]
    while len(Znuc) < NZnuc:
      Znuc_new = filter(None, data[line_znuc].split(' '))
      Znuc_new = [float(z) for z in Znuc_new]
      Znuc.extend(Znuc_new)
      line_znuc = line_znuc + 1
      
    build = []
    for i in range(N):
      Z = [Znuc[Zind[i]-1]]
      Z.extend(R[i])
      build.append(Z)
    self.molecule = qtk.Molecule()
    self.molecule.build(build)

    if self.scf_step > 0:
      fStr = filter(lambda x: 'tesian forces (hartree/b' in x, data)[-1]
      fInd = data.index(fStr)
      fData = data[fInd+1:fInd+1+N]
      force = []
      for f in fData:
        fStr = filter(None, f.split(' '))[1:]
        force.append([float(fs) for fs in fStr])
      self.force = np.array(force)

    self.occupation = []
    try:
      r1p = re.compile(r'^[ a-z]{17} +[ 0-9.E+-]+$')
      r2p = re.compile(r'^ +[a-z]+ +.*$')
      report = filter(r2p.match, filter(r1p.match, data))
      occ_pattern = filter(lambda x: ' occ ' in x, report)[-1]
      occ_pattern_ind = len(report) - report[::-1].index(occ_pattern)
      occ_pattern_end = report[occ_pattern_ind]
      occ_ind_start = len(data) - data[::-1].index(occ_pattern) - 1
      occ_ind_end = len(data) - data[::-1].index(occ_pattern_end) - 1
      for i in range(occ_ind_start, occ_ind_end):
        for occ in filter(None, data[i].split(' ')):
          try:
            self.occupation.append(float(occ))
          except Exception as err:
            pass
    except Exception as err:
      qtk.warning("error when extracting occupation number with" +\
        " error message: %s" % str(err))

    cell_pattern = re.compile(r'^ R.*=.* G.*=.*$')
    cell_list = filter(cell_pattern.match, data)
    cell = []
    for cstr in cell_list:
      cell.append(
        [float(c) for c in filter(None, cstr.split(' '))[1:4]])
    self.lattice = np.array(cell) / 1.889726124993
    self.celldm = qtk.lattice2celldm(self.lattice)
    self.molecule.R_scale = qtk.xyz2fractional(
      self.molecule.R, self.celldm)

    eigStr = os.path.join(os.path.split(qmout)[0], '*_EIG')
    eigFileList = sorted(glob.glob(eigStr))
    if len(eigFileList) != 0:
      if 'eig_index' in kwargs:
        eig_i = kwargs['eig_index']
      else:
        eig_i = -1
      if len(eigFileList) > 1:
        qtk.warning(
          "more than one o_EIG files found " + \
          "loading last file with name: %s" % eigFileList[eig_i]
        )
      eigFile = open(eigFileList[eig_i])
      eigData = eigFile.readlines()
      eigFile.close()
      unitStr = filter(lambda x: 'Eigenvalues' in x, eigData)[0]
      unitStr = unitStr.replace('(', '').replace(')', '')
      unit = filter(None, unitStr.split(' '))[1]
      spinList = filter(lambda x: 'SPIN' in x, eigData)
      if len(spinList) != 0:
        spinFactor = 2
        maxInd = eigData.index(spinList[-1])
      else:
        spinFactor = 1
        maxInd = len(eigData)
      ind = []
      for kStr in filter(lambda x: 'kpt#' in x, eigData):
        ind.append(eigData.index(kStr))
      band = []
      kpoints = []
      if spinFactor == 1:
        for i in range(len(ind)):
          wcoord = eigData[ind[i]].split('wtk=')[-1].split(', kpt=')
          weight = float(wcoord[0])
          cStr = filter(None, wcoord[1].split('(')[0].split(' '))
          coord = [float(c) for c in cStr]
          coord.append(weight)
          kpoints.append(coord)
          s = ind[i] + 1
          if i < len(ind) - 1:
            e = ind[i+1]
          else:
            e = len(eigData)
          eig_i = filter(None, ''.join(eigData[s:e]).split(' '))
          band.append([qtk.convE(float(ew), '%s-eV' % unit)[0]
                       for ew in eig_i])
  
        self.band = np.array(band)
        self.kpoints = np.array(kpoints)
        self.mo_eigenvalues = np.array(band[0]).copy()
        if len(self.occupation) > 0:
          diff = np.diff(self.occupation)
          ind = np.array(range(len(diff)))
          pos = diff[np.where(abs(diff) > 0.5)]
          mask = np.in1d(diff, pos)
          if len(ind[mask]) > 0:
            N_state = ind[mask][0]
            vb = max(self.band[:, N_state])
            cb = min(self.band[:, N_state + 1])
            vb_pos = np.argmax(self.band[:, N_state])
            cb_pos = np.argmin(self.band[:, N_state + 1])
            self.Eg = cb - vb
            if vb_pos == cb_pos:
              self.Eg_direct = True
            else:
              self.Eg_direct = False
            self.fermi_index = N_state
  
      else:
        qtk.warning("spin polarized band data " +\
                    "extraction is not yet implemented")
    else:
      qtk.warning('no k-point information (o_EIG file) found')

  def unfold(self, k_path, folds, epsilon = 1E-5, WFK=None, overwrite=False, shift = None):

    if not WFK:
      path = self.path
      if not path:
        path = '.'
      WFK_card = '%s/*_WFK' % path
      WFK_list = sorted(glob.glob(WFK_card))
      if len(WFK_list) > 0:
        WFK = WFK_list[-1]
      else:
        qtk.exit('wavefunction file not found.')
      
    path, name = os.path.split(WFK)
    if not path:
      path = '.'
    cwd = os.getcwd()
    os.chdir(path)

    f2b_list = sorted(glob.glob('*.f2b'))
    if not f2b_list or overwrite:
      exe = qtk.setting.abinit_f2b_exe
      log_name = '%s_f2b.log' % self.name
      log = open(log_name, 'w')
      fold_str = [str(f) for f in folds]
      cmd_str = "%s %s %s" % (exe, name, ':'.join(fold_str))
      run = sp.Popen(cmd_str, shell=True, stdout=log)
      run.wait()
      log.close()
      f2b_list = sorted(glob.glob('*.f2b'))

    f2b = f2b_list[0]

    k_path = np.asarray(k_path)

    def coordTransform(V, G):
      W = np.zeros(V.shape)
      GT = G.T
      for i in range(len(V)):
          W[i,:] = GT[0,:]*V[i,0] + GT[1,:]*V[i,1] + GT[2,:]*V[i,2]
      return W

    def dp2l(X0,X1,X2):
      eps = 0.001
      denom = X2 - X1
      denomabs = np.linalg.norm(denom)
      if denomabs < eps:
          return False
      numer = np.cross( X0-X1 , X0-X2 )
      numerabs = np.linalg.norm(numer)
      return numerabs / denomabs

    _data = np.loadtxt(f2b)
    os.chdir(cwd)
    KEIG = _data[:, :3]
    EIG = np.array(_data[:, 3]) * qtk.convE(1, 'Ha-ev')[0]
    if not shift:
      EIG = EIG - np.max(self.band[:, self.fermi_index])
    else:
      EIG = EIG - shift
    EIG = EIG.tolist()
    W = _data[:, 4]

    L = []
    ENE = []
    WGHT = []
    G = self.G_lattice.copy()

    for i in range(3):
      G[i,:] = G[i,:] * folds[i]

    k_path = coordTransform(k_path, G)
    KEIG = coordTransform(KEIG, G)

    dl = 0
    for ikp in range(len(k_path) - 1):
      itr = 0
      for j in range(len(EIG)):
        itr += 1
        dist = dp2l(KEIG[j,:], k_path[ikp, :], k_path[ikp+1, :])
        B = k_path[ikp,:] - k_path[ikp+1, :]
        dk = np.linalg.norm(B)
        if dist < epsilon:
          A = k_path[ikp, :] - KEIG[j,:]
          x = np.dot(A, B) / dk
          if x > 0 and x-dk < epsilon:
            L.append(x+dl)
            ENE.append(EIG[j])
            WGHT.append(W[j])
      dl = dl + dk
    L = np.array(L)
    ENE = np.array(ENE)
    WGHT = np.array(WGHT)

    return L, ENE, WGHT
