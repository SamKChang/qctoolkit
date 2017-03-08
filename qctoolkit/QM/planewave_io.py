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
from itertools import permutations

import pkgutil
ase_eggs_loader = pkgutil.find_loader('ase')
ase_found = ase_eggs_loader is not None
if ase_found:
  import ase.dft.kpoints as asekpt
spg_eggs_loader = pkgutil.find_loader('spglib')
spg_found = spg_eggs_loader is not None
if spg_found:
  import spglib

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

  def vbEdge(self):

    if hasattr(self, 'band_symmetrized'):
      band = self.band_symmetrized
    else:
      band = self.band

    diff = np.diff(self.occupation)
    pos = diff[np.where(abs(diff) > 0.5)]
    mask = np.in1d(diff, pos)
    ind = np.array(range(len(diff)))
    N_state = ind[mask][0]
    vb = max(band[:, N_state])
    return vb, N_state

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

  def _fillBox(self, krange, n_point):
    if not hasattr(self, '_filled_band'):
      lattice = self._k_basis()
  
      kpoints = copy.deepcopy(self.kpoints)
      band = copy.deepcopy(self.band)
  
      for i in [-1, 0, 1]:
        for j in [-1, 0, 1]:
          for k in [-1, 0, 1]:
            if i**2 + j**2 + k**2 != 0:
              vec = lattice[0]*i + lattice[1]*j + lattice[2]*k
              kpoints = np.vstack([kpoints, self.kpoints + 2*vec])
              band = np.vstack([band, self.band])
      self._filled_band = band
      self._filled_kpoints = kpoints

      self._kmeshgrid = np.mgrid[
        -krange:krange:n_point,
        -krange:krange:n_point,
        -krange:krange:n_point]
      self._rmeshgrid = [
        np.mgrid[-krange:krange:n_point],
        np.mgrid[-krange:krange:n_point],
        np.mgrid[-krange:krange:n_point],
      ]

  def _interpolate(self, single_band, n_point=13j, krange=1):


    self._fillBox(krange, n_point)
    band = self._filled_band
    kpoints = self._filled_kpoints
    
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
      tuple(self._kmeshgrid), 
      method='linear'
    )

    interpolated = RegularGridInterpolator(
      tuple(self._rmeshgrid), inter_k - vb
    )
    return interpolated, cb

  def _kPath(self, special_points, dk):
    spk = self.special_kpoints
    out = []
    pos = 0
    tick_pos = [0]
    tick_txt = []
    for i in range(len(special_points)-1):
      ci = special_points[i]
      cf = special_points[i+1]
      if ci not in spk:
        qtk.exit("special point %c not in reconized" % ci)
      if cf not in spk:
        qtk.exit("special point %c not in reconized" % cf)
      if len(tick_txt) == 0:
        tick_txt.append("$\mathrm{" + ci + "}$")
      Ri = spk[ci]
      Rf = spk[cf]
      if len(out) == 0:
        out.append(Ri)
      vector = np.array(Rf) - np.array(Ri)
      n_points = int(ceil(np.linalg.norm(vector) / float(dk)))
      pos = pos + n_points
      tick_pos.append(pos)
      tick_txt.append("$\mathrm{" + cf + "}$")
      for i in range(1, n_points+1 ):
        point = list(np.round(
          np.array(Ri) + (i/float(n_points))*vector, decimals=3
        ))
        out.append(point)
      tick_txt = [
        '$\Gamma$' if x == '$\mathrm{G}$' else x for x in tick_txt]
    return out, [tick_pos, tick_txt]

  def plot_band(self, path, dk=0.1, n_cb=2, n_vb=2, krange=1, ax=None):

    if not hasattr(self, 'kpoints_symmetrized'):
      qtk.warning('Make sure full kmesh is available. Brillouinize might help')

    p, ticks = self._kPath(path, dk)
    x = np.linspace(0, len(p)-1, 100)
    tick_pos, tick_txt = ticks[0], ticks[1]
    if ax is None:
      fig, ax = plt.subplots()
    bands = []
    bands_out = []
    occ = np.array(self.occupation)
    vb_top = np.where(occ < 1)[0][0]
    for b in range(vb_top - n_vb, vb_top + n_cb):
      bout, cb = self._interpolate(b, krange=krange)
      band_points = bout(p)
      bands_out.append(bout)
      interb = interp1d(range(len(p)), band_points, kind='linear')
      if b < vb_top:
        color = 'r'
      else:
        color = 'b'
      if b == vb_top or b == vb_top + 1:
        bands.append(bout)
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

  def brillouinScale(self, scale):
    length = np.linalg.norm(self.special_kpoints['X'])
    factor = float(scale) / length
    for k, v in self.special_kpoints.items():
      self.special_kpoints[k] = (np.array(v) * factor).tolist()

  def brillouinize(self, kmesh=None, kgrid=None):

    if not hasattr(self, 'kpoints') or len(self.kpoints) == 0:
      qtk.exit('kpoints information not available')
    else:
      if not hasattr(self, 'kpoints_symmetrized'):
        k_old = self.kpoints[:,:3]
        b_old = self.band
      else:
        k_old = self.kpoints_symmetrized[:,:3]
        b_old = self.band_symmetrized

    if kgrid is None and kmesh is not None:
      try:
        kgrid, new_band = self._spg_grid(k_old, b_old, kmesh)
      except:
        try:
          kgrid, new_band = self._ase_grid(k_old, b_old, kmesh)
        except Exception as e:
          qtk.exit('kpoint generation failed: %s' % e)
    elif kmesh is None and kgrid is not None:
      new_band = self._kgrid_template(kgrid)

    if not hasattr(self, 'kpoints_symmetrized'):
      self.kpoints_symmetrized = self.kpoints
      self.band_symmetrized = self.band
    self.kpoints = kgrid
    self.band = new_band

  def _spg_grid(self, k_old, b_old, kmesh):

    shift_ind = np.argmin(
      np.linalg.norm(np.abs(k_old[:,:3]), axis=1)
    )
    shift = k_old[shift_ind]
    
    lattice = self.lattice
    positions = self.molecule.R_scale
    numbers = self.molecule.Z
    cell = (lattice, positions, numbers)
    mapping, grid = spglib.get_ir_reciprocal_mesh(
      kmesh, cell, is_shift=shift)
    kgrid = grid.astype(float) / kmesh
    kgrid_shifted = kgrid + np.array([1,1,1])

    groups_dict = {}
    for i in range(len(mapping)):
      k_id = mapping[i]
      if k_id not in groups_dict:
        groups_dict[k_id] = [i]
      else:
        groups_dict[k_id].append(i)
    groups = list(groups_dict.values())
    
    
    new_band = np.zeros([len(kgrid), len(b_old[0])])
    
    for i in range(len(k_old[:,:3])):
      k = k_old[i,:3]
      norm = []
      for k_permute in permutations(k):
        norm.append(np.linalg.norm(kgrid - k_permute, axis=1))
        norm.append(np.linalg.norm(kgrid + k_permute, axis=1))
      for k_permute in permutations(k):
        norm.append(np.linalg.norm(kgrid - np.abs(k_permute), axis=1))
        norm.append(np.linalg.norm(kgrid + np.abs(k_permute), axis=1))
      norm = np.min(np.stack(norm), axis=0)
      try:
        key = np.where(norm < 1E-3)[0][0]
      except Exception as e:
        msg = 'kpoint miss match with error message: ' + str(e)
        msg = msg + '\ncurrtent kpoint: ' + str(k)
        qtk.exit(msg)
      for g in groups:
        if key in g:
          for key_g in g:
            new_band[key_g] = b_old[i]

    return kgrid, new_band

#
#    if kgrid is None and kmesh is not None:
#      new_kpoints, k_dup, ind_key = self._ase_grid(kmesh, k_old)
#      new_band = []
#      for k_new in new_kpoints:
#        norm1 = np.linalg.norm(k_dup - k_new, axis=1)
#        norm2 = np.linalg.norm(k_dup + k_new, axis=1)
#        try:
#          key = np.where(norm1 < 1E-6)[0][0]
#        except:
#          try:
#            key = np.where(norm2 < 1E-6)[0][0]
#          except Exception as e:
#            print k_new
#            qtk.exit('kmesh not matched with error message: %s' % str(e))
#        ind = ind_key[key]
#        new_band.append(b_old[ind])
#      
#      new_band = np.stack(new_band)
#
#    elif kgrid is not None:
#      new_kpoints = kgrid[:, :3]
#      k_dup = self.kpoints[:, :3]
#      ind_key = range(len(self.kpoints))
#
#      new_band = np.zeros([len(new_kpoints), len(self.band[0])])
#
#      for i in range(len(k_dup)):
#        k = k_dup[i]
#        b = b_old[i]
#        norm = np.linalg.norm(new_kpoints - k, axis=1)
#        try:
#          key = np.where(norm < 1E-4)[0][0]
#        except Exception as e:
#          print k
#          qtk.exit('kmesh not matched with error message: %s' % str(e))
#        new_band[key] = b
#
#    if not hasattr(self, 'kpoints_symmetrized'):
#      self.kpoints_symmetrized = self.kpoints
#      self.band_symmetrized = self.band
#    self.kpoints = new_kpoints
#    self.band = new_band
#

  def BZ_integrate(self, integrand):
    results = []
    integrand = np.asarray(integrand)
    assert len(integrand) == len(self.kpoints)
    w = np.asarray(self.kpoints[:, -1])
    if len(integrand.T[0]) == len(self.kpoints):
      for vec in integrand.T:
        results.append(w.dot(vec))
    else:
      results = w.dot(integrand)
    return results

  def DOS(self, sigma=0.2, E_list=None, raw=False):
    dos_raw = []
    for b in self.band.T:
      dos_raw.extend(zip(b, self.kpoints[:,-1]))
    dos_raw = np.array(dos_raw)

    if raw:
      return raw
    else:
      if not E_list:
        E_min, E_max = min(dos_raw[:, 0]), max(dos_raw[:, 0])
        p_max = E_max + 0.1 * (E_max - E_min)
        p_min = E_min - 0.1 * (E_max - E_min)
        E_list = np.linspace(p_min, p_max, 500)

      y = np.zeros(len(E_list))
      for entry in dos_raw:
        E, w = entry
        y = y + w * np.exp(-((E-E_list)/sigma)**2)\
          / (sigma * np.sqrt(np.pi))
      return np.array([E_list, y]).T

  def _kgrid_template(self, kgrid):
    new_kpoints = kgrid[:, :3]
    k_dup = self.kpoints[:, :3]
    ind_key = range(len(self.kpoints))

    new_band = np.zeros([len(new_kpoints), len(self.band[0])])

    for i in range(len(k_dup)):
      k = k_dup[i]
      b = b_old[i]
      norm = np.linalg.norm(new_kpoints - k, axis=1)
      try:
        key = np.where(norm < 1E-4)[0][0]
      except Exception as e:
        print k
        qtk.exit('kmesh not matched with error message: %s' % str(e))
      new_band[key] = b
    return new_band

  def _ase_grid(self, k_old, b_old, kmesh):
    if not ase_found:
      qtk.exit('ase not found')
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
    zeros = np.zeros(len(k_x))
    ind_key_base = range(len(k_x))
  
    k_tmp = []
    k_tmp.append(np.stack([ k_x,  k_y,  k_z]).T)
    k_tmp.append(np.stack([ k_x,  k_y, -k_z]).T)
    k_tmp.append(np.stack([ k_x, -k_y,  k_z]).T)
    k_tmp.append(np.stack([ k_x, -k_y, -k_z]).T)
    k_tmp.append(np.stack([-k_x,  k_y,  k_z]).T)
    k_tmp.append(np.stack([-k_x,  k_y, -k_z]).T)
    k_tmp.append(np.stack([-k_x, -k_y,  k_z]).T)
    k_tmp.append(np.stack([-k_x, -k_y, -k_z]).T)
    
    ind_key_tmp = np.concatenate([ind_key_base for _ in k_tmp])
    ind_key = ind_key_tmp
    ind_key_base = ind_key.copy()
    k_dup = np.concatenate(k_tmp)
  
    k_x = k_dup[:,0]
    k_y = k_dup[:,1]
    k_z = k_dup[:,2]
    
    k_tmp = []
    k_tmp.append(np.stack([k_x, k_y, k_z]).T)
    k_tmp.append(np.stack([k_z, k_x, k_y]).T)
    k_tmp.append(np.stack([k_y, k_z, k_x]).T)
    k_tmp.append(np.stack([k_z, k_y, k_x]).T)
    k_tmp.append(np.stack([k_x, k_z, k_y]).T)
    k_tmp.append(np.stack([k_y, k_x, k_z]).T)
    
    ind_key_tmp = np.concatenate([ind_key_base for _ in k_tmp])
    ind_key = np.concatenate([ind_key, ind_key_tmp])
    ind_key_base = ind_key.copy()
  
    k_tmp = np.concatenate(k_tmp)
    k_dup = np.vstack([k_dup, k_tmp])

    new_band = np.zeros([len(new_kpoints), len(self.band[0])])

    for i in range(len(k_dup)):
      k = k_dup[i]
      b = b_old[i]
      norm = np.linalg.norm(new_kpoints - k, axis=1)
      key = np.where(norm < 1E-3)[0][0]
      new_band[key] = b

    return new_kpoints, new_band
