import glob, copy, os
import qctoolkit as qtk
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import interp1d

class PathData(qtk.QMData):
  """
  analysis on alchemical path results
  """
  _file_list = []
  _cube_list = []

  @classmethod
  def loadCubeList(cls, path_list, program=qtk.setting.qmcode):
    if program == 'cpmd':
      cls._file_list = path_list
      _para = [[name] for name in cls._file_list]
      if len(_para)<3*qtk.setting.cpu_count:
        cls._cube_list = qtk.parallelize(qtk.CUBE, _para, 
                                     block_size=1)
      else:
        cls._cube_list = qtk.parallelize(qtk.CUBE, _para)
    else:
      qtk.exit("density of alchemical path is "\
               +"not yet implemented for %s" % self.program)
    
  @classmethod
  def loadAllCube(cls, path, pattern, 
                  program=qtk.setting.qmcode):
    """
    centralized parallel CUBE loading
    """
    if program == 'cpmd':
      cls._file_list = sorted(glob.glob(
        path + "/" + pattern + "/*.cube"))
      _para = [[name] for name in cls._file_list]
      if len(_para)<3*qtk.setting.cpu_count:
        cls._cube_list = qtk.parallelize(qtk.CUBE, _para, 
                                     block_size=1)
      else:
        cls._cube_list = qtk.parallelize(qtk.CUBE, _para)
    else:
      qtk.exit("density of alchemical path is "\
               +"not yet implemented for %s" % self.program)

  def __init__(self, path, pattern, program, **kwargs):
    qtk.QMData.__init__(self, path, pattern, program, **kwargs)
    self.cube_list = []

  def loadCube(self):
    # cpmd implemetation
    if self.program == 'cpmd':
      self.cube_name = sorted(glob.glob(
        self.path + "/" + self.pattern + "/*.cube"))
      if len(PathData._cube_list) > 0:
        _clist = PathData._cube_list
        _nlist = PathData._file_list
        self.cube_list = [_clist[_nlist.index(name)]\
                          for name in self.cube_name]
      else:
        _para = [[name] for name in self.cube_name]
        self.cube_list = qtk.parallelize(qtk.CUBE, _para, 
                                           block_size=1)
    else:
      qtk.exit("density of alchemical path is "\
               +"not yet implemented for %s" % self.program)

  def shiftCube(self, vector):
    for cube in self.cube_list:
      cube.shift(vector)

  def plotCube(self, ax=None, **kwargs):
    if ax is None:
      if 'name' in kwargs:
        _name = kwargs['name']
      else:
        _name = 'density_plot'
      fig = plt.figure(_name)
      ax = fig.gca(projection='3d')
      if 'no_show' in kwargs and kwargs['no_show']:
        _show = False
      else:
        _show = True
    else:
      _show = False
    if 'spacing' in kwargs:
      spacing = kwargs['spacing']
    else:
      spacing = 1
    if 'axis' in kwargs:
      _axis = kwargs['axis']
    else:
      _axis = 0
    if 'track' in kwargs:
      track = kwargs['track']
    else:
      track = []
      
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)

    # get values from first cube file
    x, z = self.cube_list[0].plot(no_show=True, axis=_axis)
    if 'plot_range' in kwargs:
      _min = kwargs['plot_range'][0]
      _max = kwargs['plot_range'][1]
    else:
      _min = x[0]
      _max = x[-1]
    l = np.linspace(0, 1, len(self.cube_list))
    X, Y = np.meshgrid(x, l)
    Z = np.atleast_2d(z)
    # load every cube file to meshgrid for contour plot
    for ind in range(1, len(self.cube_list)):
      x, z = self.cube_list[ind].plot(no_show=True, axis=_axis)
      Z = np.vstack([Z, z])

    itr = 0
    cube_track = []
    for l in np.linspace(0, 1, len(self.cube_list)):
      x, z = self.cube_list[itr].plot(no_show=True, axis=_axis)
      y = [l for _ in range(len(x))]
      # contruct tracking data points
      if len(track) > 0:
        cube_interp = interp1d(x, z)
        cube_track.append(cube_interp(track))
      # only selected cubes are plotted
      if itr % spacing == 0:
        # chop off extra data
        x[x < _min] = np.nan
        x[x > _max] = np.nan
        # plot integrated density according to specified axis
        ax.plot(x, y, z, color = 'blue')
      itr += 1
    cube_track = np.array(cube_track).T

    # construct tracking coordinates
    if len(track) > 0:
      if 'track_style' in kwargs:
        track_style = kwargs['track_style']
      else:
        track_style = ['--' for _ in track]
      track_y = np.linspace(0, 1, len(self.cube_list))
      for i in range(len(cube_track)):
        track_z = cube_track[i]
        track_x = np.array([track[i] for x in range(len(track_z))])
        # plot tracking lines
        ax.plot(track_x, track_y, track_z, 
                color='red', linestyle=track_style[i])
        if 'track_project' in kwargs and kwargs['track_project']:
          proj_x = np.array([_min for x in range(len(track_z))])
          ax.plot(proj_x, track_y, track_z, 
                  color='red', linestyle=track_style[i])
       

    if 'levels' in kwargs:
      levels = kwargs['levels']
    else:
      lmin = np.min(Z) - (np.max(Z) - np.min(Z))/20
      lmax = np.max(Z) + (np.max(Z) - np.min(Z))/20
      #levels = np.logspace(0.05, np.log10(np.max(Z)), num=10) - 1
      levels = np.linspace(lmin, lmax, 10)
    if 'contour_floor' in kwargs:
      cfloor = kwargs['contour_floor']
    else:
      cfloor = kwargs['contour_floor']
    floor = np.min(Z) - (np.max(Z) - np.min(Z))/20
    ceil = np.max(Z) + (np.max(Z) - np.min(Z))/20
    cset = ax.contour(X, Y, Z, levels,
                      zdir='z', offset=cfloor, colors='black')

    if 'xlabel' in kwargs:
      xlabel = kwargs['xlable']
    else:
      xlabel = "$x$ [$\AA$]"
    if 'ylabel' in kwargs:
      ylabel = kwargs['ylable']
    else:
      ylabel = "$\lambda$"
    if 'zlabel' in kwargs:
      zlabel = kwargs['zlable']
    else:
      zlabel = r"$\rho(x, \lambda)\ [a.u.]$"

    ax.set_xlim3d(_min, _max)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylim3d(0,1)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_zlim3d(floor, ceil)
    ax.set_zlabel(zlabel, 
                  fontsize=16, rotation=90)

    if 'legend' in kwargs:
      if 'legend_size' in kwargs:
        legend_size = kwargs['legend_size']
      else:
        legend_size = 12
      ax.text2D(0.7, 0.9, kwargs['legend'],
                fontsize=legend_size,
                transform=ax.transAxes,
                bbox={'facecolor':'white', 'alpha':1, 'pad':10})
    if _show:
      plt.show()

  # add cube list difference to __add__
  def __add__(self, other):
    _out = super(PathData, self).__add__(other)
    _out.cube_list = []
    if len(self.cube_list) == len(other.cube_list):
      for i in range(len(self.cube_list)):
        _out.cube_list.append(\
          self.cube_list[i] + other.cube_list[i])
    return _out

  # add cube list difference to __sub__
  def __sub__(self, other):
    _out = super(PathData, self).__sub__(other)
    _out.cube_list = []
    if len(self.cube_list) == len(other.cube_list):
      for i in range(len(self.cube_list)):
        _out.cube_list.append(\
          self.cube_list[i] - other.cube_list[i])
    return _out

class PathScan(object):
  """
  generate input files to scan alchemical path
  """
  def __init__(self, xyz_i, xyz_f, 
               program=qtk.setting.qmcode, **kwargs):
    if 'name' in kwargs:
      self.name = kwargs['name']
    elif type(xyz_i)==str and type(xyz_f)==str:
      ref_stem, _ = os.path.splitext(xyz_i)
      tar_stem, _ = os.path.splitext(xyz_f)
      self.name = ref_stem + '-' + tar_stem
    else:
      self.name = 'A2B'
    self.program = program

    if type(xyz_i) == str:
      self.initial = qtk.Molecule()
      self.initial.read(xyz_i)
    elif type(xyz_i) == qtk.Molecule:
      self.initial = xyz_i
    if type(xyz_f) == str:
      self.final = qtk.Molecule()
      self.final.read(xyz_f)
    elif type(xyz_f) == qtk.Molecule:
      self.final = xyz_f

    self.inp_base = qtk.QMInp(xyz_i, program)
    if 'inp_base' in kwargs:
      if kwargs['inp_base'] is not None:
        inp_tmp = copy.deepcopy(kwargs['inp_base'])
        self.inp_base.inp.setting = inp_tmp.inp.setting

    # !!!!!!!!!!!!!!!!!!!!!!!
    # construct mutation list
    # tar_list: list of matching location
    tar_list = [i for i in range(self.final.N)]
    self.mutation_ind = []
    self.mutation_ref = []
    self.mutation_tar = []
    self.mutation_crd = []
    # loop through every atom in initial system
    for i in range(self.initial.N):
      Zi = self.initial.type_list[i]
      Ri = self.initial.R[i]
      to_void = True
      # loop through every atom in final system
      for j in range(self.final.N):
        # every atom-pair distance
        Rj = self.final.R[j]
        diff = np.linalg.norm(Ri-Rj)
        # search for atoms at same location
        if diff < 10E-6:
          to_void = False
          tar_list[j] = False
          Zj = self.final.type_list[j]
          if Zi != Zj:
            self.mutation_ind.append(i)
            self.mutation_ref.append(Zi)
            self.mutation_tar.append(Zj)
            self.mutation_crd.append(Ri)
      # when no distance matched, set atom target to void
      if to_void:
        self.mutation_ind.append(i)
        self.mutation_ref.append(Zi)
        self.mutation_tar.append('VOID')
        self.mutation_crd.append(Ri)

    # loop through every target atom for non-matching atoms
    N = self.initial.N
    for j in range(self.final.N):
      if tar_list[j]:
        self.mutation_ind.append(N)
        self.mutation_ref.append('VOID')
        self.mutation_tar.append(self.final.type_list[j])
        self.mutation_crd.append(self.final.R[j])

  def PPString(self, ind):
    ref = self.mutation_ref[ind]
    tar = self.mutation_tar[ind]
    if ref == 'VOID':
      ref = 'V'
    if tar == 'VOID':
      tar = 'V'
    return ref+str(2)+tar 

  def l(self, l_in):
    if self.program == 'cpmd':
      ppstr = "_%03d.psp" % (l_in*100)
      out = copy.deepcopy(self.inp_base)
      for i in range(len(self.mutation_ind)):
        ind = self.mutation_ind[i]
        pp = self.PPString(i)+ppstr
        if ind < out.inp.structure.N:
          out.setAtom(ind+1, pp)
        else:
          out.addAtom(pp, self.mutation_crd[i])
      return out
    else:
      qtk.exit("program %s is not implemented for PathScan"\
               % self.program)
    
  def fullPath(self, stem, dl=0.1):
    inp_list = []
    for l in np.linspace(0,1,(1/dl)+1):
      inp_list.append(self.l(l))
    return inp_list
