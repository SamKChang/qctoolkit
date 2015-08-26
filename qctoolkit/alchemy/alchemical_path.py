import glob
import qctoolkit as qtk
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import interp1d

class PathData(qtk.QMData):
  def __init__(self, path, pattern, program, **kwargs):
    qtk.QMData.__init__(self, path, pattern, program, **kwargs)

  def loadCube(self):
    # cpmd implemetation
    if self.program == 'cpmd':
      self.cube_name = sorted(glob.glob(
        self.path + "/" + self.pattern + "/*.cube"))
      _para = [[name] for name in self.cube_name]
      self.cube_list = qtk.parallelize(qtk.CUBE, _para, 
                                       block_size=1)
    else:
      qtk.exit("density of alchemical path is "\
               +"not yet implemented for %s" % self.program)

  def plotCube(self, ax=None, **kwargs):
    if ax is None:
      fig = plt.figure('density_plot')
      ax = fig.gca(projection='3d')
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
    x, z = self.cube_list[0].plot(no_show=True)
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
      track_y = np.linspace(0, 1, len(self.cube_list))
      for i in range(len(cube_track)):
        track_z = cube_track[i]
        track_x = np.array([track[i] for x in range(len(track_z))])
        # plot tracking lines
        ax.plot(track_x, track_y, track_z, 
                color='red', linestyle='--')

    if 'levels' in kwargs:
      levels = kwargs['levels']
    else:
      levels = np.logspace(0.05, np.log10(np.max(Z)), num=10) - 1
    floor = np.min(Z) - (np.max(Z) - np.min(Z))/20
    ceil = np.max(Z) + (np.max(Z) - np.min(Z))/20
    cset = ax.contour(X, Y, Z, levels,
                      zdir='z', offset=floor, colors='black')

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
      ax.text2D(0.7, 0.9, kwargs['legend'],
                transform=ax.transAxes,
                bbox={'facecolor':'white', 'alpha':1, 'pad':10})
