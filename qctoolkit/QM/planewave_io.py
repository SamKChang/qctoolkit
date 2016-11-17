import qctoolkit as qtk
from general_io import GenericQMInput
from general_io import GenericQMOutput
import universal as univ
import numpy as np
import copy
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
from math import ceil
import matplotlib.pyplot as plt

import pkgutil
ase_eggs_loader = pkgutil.find_loader('ase')
ase_found = ase_eggs_loader is not None
if ase_found:
  import ase.dft.kpoints as asekpt

class PlanewaveInput(GenericQMInput):
  """
  From PlanwaveInput:
  generic class holder for plane wave qmcode. It provide basic
  default settings.
  """
  __doc__ = GenericQMInput.__doc__ + __doc__
  def __init__(self, molecule, **kwargs):
    GenericQMInput.__init__(self, molecule, **kwargs)

    self.setting.update(kwargs)

    if 'cutoff' not in kwargs:
      self.setting['cutoff'] = 100
    if not self.setting['periodic'] and 'isolation' not in kwargs:
      self.setting['isolation'] = 'mt'
    self.pp_files = []
    if 'periodic' in self.setting and self.setting['periodic']:
      self.celldm2lattice()
    if 'pp_type' not in kwargs:
      self.setting['pp_type'] = 'geodecker'
    if 'full_kmesh' not in self.setting:
      self.setting['full_kmesh'] = False
    if 'theory' in self.setting and self.setting['theory'] == 'hse06':
      if 'pp_theory' not in self.setting:
        self.setting['pp_theory'] = 'pbe'
    if 'fractional_coordinate' not in kwargs:
      self.setting['fractional_coordinate'] = False

    univ.getCelldm(self) 

  def write(self, name=None, **kwargs):
    if self.setting['periodic']:
      self.celldm2lattice()
    inp, molecule = \
      GenericQMInput.write(self, name, **self.setting)

    if 'pp_list' in self.setting:
      pp_list = self.setting['pp_list']
      itr = 1
      for pp_data in pp_list:
        pp_inds = pp_data[0]
        if type(pp_inds) is not list:
          pp_inds = [pp_inds]
        pp_name = pp_data[1]
        pp = pp_data[2]
        pp.setting['program'] = self.setting['program']
        pp.write(pp_name, inplace=False)
        molecule.setAtoms(pp_inds, string=pp_name)
        Zn = molecule.type_list[pp_inds[0]]
        molecule.setAtoms(pp_inds, element=Zn + str(itr))
        itr += 1

    return inp, molecule

  def celldm2lattice(self):
    cd = self.setting['celldm']
    if 'scale' in self.setting:
      sc = self.setting['scale']
    else:
      sc = [1.0 for i in range(3)]
    self.setting['lattice'] = qtk.celldm2lattice(cd, scale=sc)

  def cornerMargin(self, *args, **kwargs):
    pass

class PlanewaveOutput(GenericQMOutput):
  def __init__(self, output=None, **kwargs):
    GenericQMOutput.__init__(self, output, **kwargs)
    self.special_kpoints = {
      'G':[0,0,0],
      'X':[1,0,0],
      'W':[1, 0.25, 0],
      'K':[0.75, 0.75, 0],
      'L':[0.5, 0.5, 0.5],
      'U':[0.1, 0.25, 0.25],
    }

  def _k_basis(self):

    basis = []
    for axis in range(3):

      i = (axis + 1) % 3
      j = (axis + 2) % 3
  
      k = self.kpoints
      k_i = k[np.where(k[:,i] == 0)]
      k_ij = k_i[np.where(k_i[:, j] == 0)]
      ind = np.argmax(k_ij[:,axis])
      basis.append(k_ij[ind])

    return np.stack(basis)


  def _interpolate(self, single_band, n_point=13j, krange=1):

    if not hasattr(self, 'kpoints_symmetrized'):
      qtk.exit('brillouinize is necessary')

    lattice = self._k_basis()

    kpoints = copy.deepcopy(self.kpoints)
    band = copy.deepcopy(self.band)
    for i in [-1, 0, 1]:
      for j in [-1, 0, 1]:
        for k in [-1, 0, 1]:
          if i**2 + j**2 + k**2 != 0:
            vec = lattice[0]*i + lattice[1]*j + lattice[2]*k
            kpoints = np.vstack([kpoints, self.kpoints + vec])
            band = np.vstack([band, self.band])

    kx,ky,kz = np.mgrid[
      -krange:krange:n_point,
      -krange:krange:n_point,
      -krange:krange:n_point]
    x = np.mgrid[-krange:krange:n_point]
    y = np.mgrid[-krange:krange:n_point]
    z = np.mgrid[-krange:krange:n_point]
    #kx, ky, kz = np.meshgrid(k,k,k)
    
    diff = np.diff(self.occupation)
    pos = diff[np.where(abs(diff) > 0.5)]
    mask = np.in1d(diff, pos)
    ind = np.array(range(len(diff)))
    N_state = ind[mask][0]
    vb = max(band[:, N_state])
    cb = min(band[:, N_state + 1])-vb

    inter_k = griddata(
      kpoints[:,:3], 
      band[:, single_band], 
      (kx, ky, kz), 
      method='linear'
    )

    return RegularGridInterpolator((x,y,z), inter_k - vb), cb

  def _kPath(self, special_points, dk = 0.01):
    spk = self.special_kpoints
    out = []
    pos = 0
    tick_pos = [0]
    tick_txt = []
    for i in range(len(special_points)-1):
      ci = special_points[i]
      cf = special_points[i+1]
      if ci not in spk:
        print "%c not in reconized" % ci
        return None
      if cf not in spk:
        print "%c not in reconized" % cf
        return None
      if len(tick_txt) == 0:
        tick_txt.append("$\mathrm{" + ci + "}$")
      Ri = spk[ci]
      Rf = spk[cf]
      if len(out) == 0:
        out.append(Ri)
      vector = np.array(Rf) - np.array(Ri)
      n_points = int(ceil(np.linalg.norm(vector) / float(dk)))
      pos = pos + n_points
      #print ci + '-' + cf + ": "  + str(n_points) + ', ' + str(pos)
      tick_pos.append(pos)
      tick_txt.append("$\mathrm{" + cf + "}$")
      for i in range(1, n_points+1 ):
        point = list(np.round(
          np.array(Ri) + (i/float(n_points))*vector, decimals=3
        ))
        #print point
        out.append(point)
      tick_txt = [
        '$\Gamma$' if x == '$\mathrm{G}$' else x for x in tick_txt]
    return out, [tick_pos, tick_txt]

  def plot_band(self, path, dk=0.1, n_cb=2, n_vb=2, krange=1, ax=None):
    p, ticks = self._kPath(path, dk)
    x = np.linspace(0, len(p)-1, 100)
    tick_pos, tick_txt = ticks[0], ticks[1]
    if ax is None:
      fig, ax = plt.subplots()
    bands = []
    occ = np.array(self.occupation)
    vb_top = np.where(occ < 1)[0][0]
    for b in range(vb_top - n_vb, vb_top + n_cb):
      bout = self._interpolate(b, krange=krange)
      band_points, cb = bout[0](p), bout[1]
      interb = interp1d(range(len(p)), band_points, kind='linear')
      if b < vb_top:
        color = 'r'
      else:
        color = 'b'
      if b == vb_top or b == vb_top + 1:
        bands.append(bout[0])
      plt.plot(x,interb(x), color=color)
      plt.plot(band_points, ls='', marker='x', color=color)
    for pos in tick_pos:
      plt.axvline(pos, color='k')
    plt.axhline(0, color='r', ls='--')
    plt.axhline(cb, color='b', ls='--')
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_txt)
    ax.set_xlim([0, len(p)-1])

    return ax

  def brillouinize(self, kmesh):
    if not ase_found:
      qtk.exit('ase not found')
    if not hasattr(self, 'kpoints') or len(self.kpoints) == 0:
      qtk.exit('kpoints information not available')
    else:
      if not hasattr(self, 'kpoints_symmetrized'):
        k_old = self.kpoints[:,:3]
        b_old = self.band
      else:
        k_old = self.kpoints_symmetrized[:,:3]
        b_old = self.band_symmetrized
    k = asekpt.monkhorst_pack(kmesh)
    k_min_ind = np.argmin(np.linalg.norm(k, axis=1))
    k_min = k[k_min_ind]
    k_shifted = k - k_min

    k_min_old_ind = np.argmin(np.linalg.norm(k_old, axis=1))
    k_min_old = k_old[k_min_old_ind]
    kpoints = k_old - k_min_old
    new_kpoints = k_shifted + k_min_old
    k_x = kpoints[:, 0]
    k_y = kpoints[:, 1]
    k_z = kpoints[:, 2]
    ind_key_base = range(len(k_x))
    ind_key = range(len(k_x))
    
    k_dup = []
    
    k_dup.append(np.stack([ k_x,  k_y,  k_z]).T)
    k_dup.append(np.stack([ k_x,  k_y, -k_z]).T)
    k_dup.append(np.stack([ k_x, -k_y,  k_z]).T)
    k_dup.append(np.stack([ k_x, -k_y, -k_z]).T)
    k_dup.append(np.stack([-k_x,  k_y,  k_z]).T)
    k_dup.append(np.stack([-k_x,  k_y, -k_z]).T)
    k_dup.append(np.stack([-k_x, -k_y,  k_z]).T)
    k_dup.append(np.stack([-k_x, -k_y, -k_z]).T)
    
    for i in range(len(k_dup)-1):
      ind_key.extend(ind_key_base)
    ind_key_base = copy.deepcopy(ind_key)
    k_dup = np.concatenate(k_dup)

    k_x = k_dup[:,0]
    k_y = k_dup[:,1]
    k_z = k_dup[:,2]
    
    k_dup = []
    k_dup.append(np.stack([k_x, k_y, k_z]).T)
    k_dup.append(np.stack([k_z, k_x, k_y]).T)
    k_dup.append(np.stack([k_y, k_z, k_x]).T)
    k_dup.append(np.stack([k_z, k_y, k_x]).T)
    k_dup.append(np.stack([k_x, k_z, k_y]).T)
    k_dup.append(np.stack([k_y, k_x, k_z]).T)
    
    for i in range(len(k_dup)-1):
      ind_key.extend(ind_key_base)
    ind_key_base = ind_key
    k_dup = np.concatenate(k_dup)
    
    new_band = []
    for k_new in new_kpoints:
      norm = np.linalg.norm(k_dup - k_new, axis=1)
      try:
        key = np.where(norm < 1E-4)[0][0]
      except Exception as e:
        qtk.exit('kmesh not matched with error message: %s' % str(e))
      ind = ind_key[key]
      new_band.append(b_old[ind])
    
    new_band = np.stack(new_band)

    if not hasattr(self, 'kpoints_symmetrized'):
      self.kpoints_symmetrized = self.kpoints
      self.band_symmetrized = self.band
    self.kpoints = new_kpoints
    self.band = new_band
