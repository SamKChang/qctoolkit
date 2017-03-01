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

      ###########
      # dataset #
      ###########
      if 'band_scan' in self.setting:
        self.content['dataset'] = odict()
        self.content['dataset']['ndtset'] = 2
        if molecule.symmetry and molecule.symmetry.lower() == 'fcc':
          if 'kshift' not in self.setting:
            self.setting['kshift'] = [
              [0.5, 0.5, 0.5],
              [0.5, 0.0, 0.0],
              [0.0, 0.5, 0.0],
              [0.0, 0.0, 0.5],
            ]

      ###################
      # restart section #
      ###################
      if 'restart' in self.setting and self.setting['restart']:
        self.content['restart'] = odict([
          ('irdwfk', 1),
          ('getwfk', -1),
          #('iscf', -2), # for non-scf calculations
        ])
      if 'restart_density' in self.setting \
      and self.setting['restart_density']:
        self.content['restart']['irdden'] = 1
        self.content['restart']['getden'] = -1
  
      ##################
      # system section #
      ##################
      if 'kmesh' not in self.setting:
        self.setting['kmesh'] = [1,1,1]
      if 'kshift' not in self.setting:
        self.setting['kshift'] = [0.0,0.0,0.0]

      self.content['system'] = odict()
      if self.setting['full_kmesh']:
        self.content['system']['kptopt'] = 3
      self.content['system']['ngkpt'] = self.setting['kmesh']
      if len(np.array(self.setting['kshift']).shape) > 1:
        self.content['system']['nshiftk'] = len(self.setting['kshift'])
      self.content['system']['shiftk'] = self.setting['kshift']
      nbnd = None
      if 'ks_states' in self.setting and self.setting['ks_states']:
        vs = int(round(molecule.getValenceElectrons() / 2.0))
        nbnd = self.setting['ks_states'] + vs
        if 'd_shell' in self.setting:
          for a in molecule.type_list:
            if a in self.setting['d_shell'] and qtk.n2ve(a) < 10:
              nbnd += 5
        self.content['system']['nband'] = nbnd
      if molecule.charge != 0:
        self.content['system']['charge='] = molecule.charge
      if molecule.multiplicity != 1:
        self.content['system']['nsppol'] = 2
        self.content['system']['occopt'] = 7
      if 'save_restart' not in self.setting \
      or not self.setting['save_restart']:
        self.content['system']['prtwf'] = 0
      if 'wf_convergence' in self.setting:
        self.content['system']['tolwfr'] = \
        self.setting['wf_convergence']

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
      if 'band_scan' in self.setting:
        bnds = self.setting['band_scan']
        if len(bnds) != 2 \
        or is_strnum(bnds[0]) \
        or len(bnds[0]) != len(bnds[1]) - 1:
          qtk.exit('band_scan format: [lst_div, lst_coord]')
        bnd_content = odict([
          ('iscf', -2),
          ('getden', -1),
          ('kptopt', -3),
          ('tolwfr', self.setting['wf_convergence']),
          ('enunit', 1),
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
  
        # compute band structrue from scratch
        if 'restart' not in self.content \
        or len(self.content['restart']) == 0\
        or 'getden' not in self.content['restart']:
          self.content['system']['prtden'] = 1
          if 'nband' in self.content['system']:
            nbnd = self.content['system'].pop('nband')
          if 'prtwf' in self.content['system']:
            self.content['system'].pop('prtwf')
          for key in self.content['system'].iterkeys():
            if key[-1] != '1':
              self.content['system'][key + '1'] = \
              self.content['system'].pop(key)
          for key in bnd_content.iterkeys():
            self.content['system'][key + '2'] = bnd_content[key]

        # compute band structure by loading wavefunction
        else:
          del self.content['dataset']['ndtset']
          keep_lst = [
            'nband',
            'tolwfr',
            'nstep',
            'ecut',
          ]
          for key in self.content['system'].iterkeys():
            if key not in keep_lst:
              del self.content['system'][key]
          for key in bnd_content.iterkeys():
            if key not in self.content['restart'].keys():
              self.content['system'][key] = bnd_content[key]
        
      ################
      # cell section #
      ################
      self.content['cell'] = odict([
        ('ecut', float(self.setting['cutoff']/2.)),
        ('nstep', self.setting['scf_step']),
      ])

      #self.content['cell'] = odict()
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
        elif  molecule.symmetry.lower() == 'fcc':
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
        ('xangst', molecule.R),
      ])

      #########################
      # write content to file #
      #########################
      inp.write('# abinit input generated by qctoolkit\n')
      for section_key in self.content.iterkeys():
        if len(self.content[section_key]) > 0:
          inp.write('\n# %s section\n' % section_key)
        for key, value in self.content[section_key].iteritems():
          if is_strnum(value):
            inp.write('%s %s\n' % (key, str(value)))
          else:
            inp.write('%s' % key)
            head = ''.join([' ' for _ in key])
            if key in ['xangst', 'xcart', 'xred']:
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
          
      if 'restart_density_path' in self.setting:
        rstdn_path = self.setting['restart_density_path']
        if os.path.exists(rstdn_path):
          dn_lst = sorted(
            glob.glob(os.path.join(rstdn_path, '*o_*DEN'))
          )
          assert len(dn_lst) > 0
          dn_src = os.path.abspath(dn_lst[-1])
          pp_files.append([dn_src, name + 'i_DEN'])

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
    eigFileList = glob.glob(eigStr)
    if len(eigFileList) != 0:
      if len(eigFileList) > 1:
        qtk.warning(
          "more than one o_EIG files found" + \
          "loading last file with name: %s" % eigFileList[-1]
        )
      eigFile = open(eigFileList[-1])
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
