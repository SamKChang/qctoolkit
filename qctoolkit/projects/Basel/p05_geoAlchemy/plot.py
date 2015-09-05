#!/usr/bin/python
import qctoolkit as qtk
import qctoolkit.alchemy as qal
import matplotlib.pylab as plt
import sys, glob
import numpy as np

name = ['p1h1','p2h1','p3h1']
scan = [1.0, 2.0]
ref_path = []
tar_path = []
opt_path = []
density_path = []

def loadData(i, j, k):
  for ref in name:
    ref_list = []
    tar_list = []
    opt_list = []
    for tar in name:
#      if ref != tar:
      if ref==name[i] and tar==name[j]:
        scan_list = []
        for s in scan:
          if s==scan[k]:
            d_str = "_d%02d_*" % (s*10)
            tname = ref+'-'+tar+d_str
            tpath = qal.PathData('.',tname,'cpmd').kcal()
            tpath.extract('.*_l(...)\.out')
            tpath.index_div(100)
            scan_list.append(tpath)
            densities = glob.glob('./'+tname+'/DENSITY.cube')
            density_path.extend(densities)
          else:
            scan_list.append(None)
        tar_list.append(scan_list)

        rname = ref+'-'+tar+'_ref_*'
        rpath = qal.PathData('.',rname,'cpmd').kcal()
        rpath.extract('.*_l(...)\.out')
        rpath.index_div(100)
        ref_list.append(rpath)
        densities = glob.glob('./'+rname+'/DENSITY.cube')
        density_path.extend(densities)

        oname = ref+'-'+tar+'_opt_*'
        opath = qal.PathData('.',oname,'cpmd').kcal()
        opath.extract('.*_l(...)\.out')
        opath.index_div(100)
        opt_list.append(opath)
        densities = glob.glob('./'+oname+'/DENSITY.cube')
        density_path.extend(densities)
      else:
        ref_list.append(None)
        tar_list.append(None)
        opt_list.append(None)
    ref_path.append(ref_list)
    tar_path.append(tar_list)
    opt_path.append(opt_list)
  #return ref_path, tar_path, opt_path

#def loadDensity():
#  for i in range(3):
#    for j in range(3):
#      ref = ref_path[i][j]
#      tar = tar_path[i][j]
#      opt = opt_path[i][j]
#      ref.loadCube()
#      tar.loadCube()
#      opt.loadCube()
#      ref.shiftCube([-7.5, -7.5, -7.5])
#      opt.shiftCube([-7.5, -7.5, -7.5])
#      tar.shiftCube([-7.5, -7.5, -7.5])
#      diff_r[i][j] = tar - ref
#      diff_o[i][j] = tar - opt
      

def pathPlot(i, j, s,**kwargs):
  name = ['p1h1','p2h1','p3h1']
  name_str = name[i] + '-' + name[j] + "-d%02d" % (scan[s]*10)
  ref = ref_path[i][j]
  tar = tar_path[i][j][s]
  opt = opt_path[i][j]
  if ref != None and tar != None and opt != None:
    ref.loadCube()
    tar.loadCube()
    opt.loadCube()
    ref.shiftCube([-7.5, -7.5, -7.5])
    opt.shiftCube([-7.5, -7.5, -7.5])
    tar.shiftCube([-7.5, -7.5, -7.5])
    diff_r = tar - ref
    diff_o = tar - opt
    dref = ref.cube_list[0].structure.distance(1,2)
    dopt = opt.cube_list[0].structure.distance(1,2)
    dtar = tar.cube_list[0].structure.distance(1,2)

    fig = plt.figure(name_str)
    ax_den_ref = fig.add_subplot(2,2,1, projection='3d')
    ax_den_tar = fig.add_subplot(2,2,2, projection='3d')
    diff_r.plotCube(ax_den_ref, spacing=4,
      name=name_str + 'ref_eq', no_show=True,
      contour_floor = -0.6,
      legend=r'$\Delta_{\rm eq} \rho$',
      legend_size = 16,
      track=[0, dref, dtar], track_style=['--', '-.', ':'],
      track_project=True
    )
    diff_o.plotCube(ax_den_tar, spacing=4,
      name=name_str + 'opt', no_show=True,
      contour_floor = -0.6,
      legend=r'$\Delta_{\rm opt} \rho$',
      legend_size = 16,
      track=[0, dopt, dtar], track_style=['--', '-.', ':'],
      track_project=True
    )
    z_tickLabel = ["  % 1.1f" % -0.4, "  % 1.1f" % 0, "  % 1.1f" % 0.4]
    ax_den_ref.set_xticks(np.arange(0,11,10))
    ax_den_ref.set_yticks(np.arange(0,1.1,0.5))
    ax_den_ref.set_zticks(np.arange(-0.4,0.5,0.4))
    ax_den_ref.set_zticklabels(z_tickLabel)
    ax_den_ref.set_xlabel('\n'+r'     $x\ [\AA]$',linespacing=0.7)
    ax_den_ref.set_ylabel('\n' + r'    $\lambda$', linespacing=0.9)
    ax_den_ref.set_zlabel('\n' + r'$P(x,\lambda)$', linespacing=2)
    ax_den_ref.set_zlim3d(-0.6,0.6)
    ax_den_ref.dist=8

    ax_den_tar.set_xticks(np.arange(0,11,10))
    ax_den_tar.set_yticks(np.arange(0,1.1,0.5))
    ax_den_tar.set_zticks(np.arange(-0.4,0.5,0.4))
    ax_den_tar.set_zticklabels(z_tickLabel)
    ax_den_tar.set_xlabel('\n'+r'     $x\ [\AA]$',linespacing=0.7)
    ax_den_tar.set_ylabel('\n' + r'    $\lambda$', linespacing=0.9)
    ax_den_tar.set_zlabel('\n' + r'$P(x,\lambda)$', linespacing=2)
    ax_den_tar.set_zlim3d(-0.6, 0.6)
    ax_den_tar.dist=8

    ax_dE = fig.add_subplot(2,2,4)
    ax_Et = fig.add_subplot(2,2,3)
    ref.data.E.plot('line', ax_Et,
                    marker='v',color='red')
    tar.data.E.plot('line', ax_Et,
                    marker='o',color='black')
    opt.data.E.plot('line', ax_Et,
                    marker='^',color='blue')
    diff_r.data.E.plot('line', ax_dE,
                       marker='v',color='red')
    diff_o.data.E.plot('line', ax_dE,
                       marker='^',color='blue')
    ax_Et.set_yticks(np.arange(-16000, -7999, 2000))
    dtar = "%3.1f" % scan[s]
    ax_Et.set_yticklabels(["-16 k", '',"-12 k",'', "-8 k"])
    ax_Et.set_xlabel(r'$\lambda$', fontsize=16)
    ax_Et.set_ylabel(r'$E_{tot}$ [kcal/mol]', fontsize=16)
    ax_Et.legend([r'$E_{tot}(d_{\rm eq})$', 
                  r'$E_{tot}($'+ dtar +r'$\ \AA)$',
                  r'$E_{tot}(d_{\rm opt})$'],
                 loc=2, numpoints=1
                )
    #ax_Et.yaxis.tick_right()
    ax_dE.yaxis.tick_right()
    ax_dE.set_yticks(np.arange(-200,201,100))
    ax_dE.set_yticklabels([-200, '',0,'',200])
    ax_dE.set_xlabel('$\lambda$', fontsize=16)
    ax_dE.set_ylabel(r'$\Delta E$ [kcal/mol]', fontsize=16)
    ax_dE.legend([r'$\Delta E($'+ dtar +r'$\ \AA, d_{\rm eq})$', 
                  r'$\Delta E($'+ dtar +r'$\ \AA, d_{\rm opt})$'],
                 loc=3, numpoints=1
                )
  if not ('no_show' in kwargs and kwargs['no_show']):
    plt.show()
  return fig, name_str


i = int(sys.argv[1])-1
j = int(sys.argv[2])-1
s = int(sys.argv[3])-1
loadData(i,j,s)
qal.PathData.loadCubeList(density_path)
out, fig_name = pathPlot(i,j,s)
out.savefig(fig_name+'.pdf')
out.savefig(fig_name+'.png')

