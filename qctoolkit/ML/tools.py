import qctoolkit as qtk
import numpy as np
import copy
import pkgutil
skl_eggs_loader = pkgutil.find_loader('pyfftw')
skl_found = skl_eggs_loader is not None
if skl_found:
  from sklearn.cross_validation import ShuffleSplit
  from sklearn.cross_validation import cross_val_score
  from sklearn.linear_model import Ridge
  from sklearn.kernel_ridge import KernelRidge

def coulomb_matrix(mol, n = -2, size = 0, 
                   sort = True, nuclear_charges = True):
  if size == 0:
    size = mol.N
  if size < mol.N:
    qtk.exit("matrix size too small")
  positions = mol.R
  if nuclear_charges:
    charges = mol.Z
  else:
    charges = np.ones(mol.N)
  differences = positions[:, np.newaxis, :] \
              - positions[np.newaxis, :, :]
  distances = np.sqrt((differences ** 2).sum(axis=-1))
  distances[distances == 0] = np.nan # replace 0 for division
  if n != 0:
    invR = (distances ** n)
  else:
    invR = distances
  invR[np.isnan(invR)] = 0 # change 0 back for getting diagonal
  diag_mask = (invR == 0).astype(int)
  charge_mask_Zij = charges[:, np.newaxis] \
                  * charges[np.newaxis, :]
  charge_mask_2p4 = 0.5 * ((charges[:, np.newaxis] \
                            * charges[np.newaxis, :]) \
                            * diag_mask) ** 1.2
  cm = invR * charge_mask_Zij + charge_mask_2p4
  if sort:
    ind = np.argsort(cm.sum(axis=-1))
    cm = cm[:, ind][ind]
  out = np.zeros([size, size])
  out[:cm.shape[0], :cm.shape[1]] = cm
  return out

def coulomb_matrices(positions, nuclear_charges = None, 
                     n = -2, sort=True):
  """
  return 3D numpy array of sorted Coulomb matrix
  """

  if nuclear_charges is None:
    shape = positions.shape
    charges = np.ones([shape[0], shape[1]])
  else:
    charges = nuclear_charges

  differences = positions[..., :, np.newaxis, :] \
              - positions[..., np.newaxis, :, :]
  distances = np.sqrt((differences ** 2).sum(axis=-1))
  distances[distances == 0] = np.nan
  if n != 0:
    invR = (distances ** n)
  else:
    invR = distances
  invR[np.isnan(invR)] = 0 # change 0 back for getting diagonal
  diag_mask = (invR == 0).astype(int)
  charge_mask_Zij = charges[..., :, np.newaxis] \
                  * charges[..., np.newaxis, :]
  # diagonal part Z**2.4 / 2
  charge_mask_2p4 = 0.5 * ((charges[..., :, np.newaxis] \
                            * charges[..., np.newaxis, :]) \
                            * diag_mask) ** 1.2
  out = invR * charge_mask_Zij + charge_mask_2p4
  if sort:
    ind = np.argsort(out.sum(axis=-1))
    # numpy fancy indexing with broadcast
    out = out[
      np.arange(len(out))[:, np.newaxis, np.newaxis], 
      ind[:, :, np.newaxis], 
      ind[:, np.newaxis, :]
    ]
  return out

def krrScore(data, 
             n_samples_list = None,
             kernel='laplacian', 
             cv = None, 
             threads = 1, 
             alphas = [1e-11],
             gammas = [0.01, 1],
             descriptor = coulomb_matrices,
             descriptor_setting = {}
            ):
  """
  return scores in the format of [alphas, gammas, samples, cv]
  """
  if n_samples_list is None:
    n_samples_list = [
      int(len(data['E']) / 10.),
      int(len(data['E']) / 5.),
      int(len(data['E']) / 2.),
    ]
  
  dsetting = copy.deepcopy(descriptor_setting)
  if 'charges' in dsetting and dsetting['charges']:
    dsetting['charges'] = data['Z']
  matrix_list = descriptor(data['xyz'], **dsetting)
  E = data['E']
  
  if cv is None:
    cv = ShuffleSplit(len(E), 
                      n_iter=5,
                      test_size=.1, 
                      random_state=42)

  qtk.report("ML.tools.krrScores setting",
             "alphas:", alphas, "\n", 
             "gammas:", gammas, "\n",
             "n_samples_list:", n_samples_list, "\n",
             "cross_validation:", cv, "\n",
             "cv_threads:", threads, "\n",
             "final score format: [alphas, gammas, samples, cv]")
      
  all_scores = []
  for alpha in alphas:
    alpha_scores = []
    all_scores.append(alpha_scores)
    for gamma in gammas:
      gamma_scores = []
      alpha_scores.append(gamma_scores)
      kernel_ridge = KernelRidge(alpha=alpha, 
                                 gamma=gamma, 
                                 kernel=kernel)
      for n_samples in n_samples_list:
        cv_ = [(train[:n_samples], test) for train, test in cv]
        scores = cross_val_score(kernel_ridge, 
                                 matrix_list.reshape(
                                   len(matrix_list), -1
                                 ), 
                                 E, 
                                 cv=cv_, 
                                 n_jobs=threads, 
                                 scoring=r'neq_mean_absolute_error')
        gamma_scores.append(scores)
              
  # [alphas, gammas, samples, cv]
  return -np.array(all_scores)

def pack(data_list, **kwargs):
  if isinstance(data_list[0], qtk.Molecule):
    typ = 'molecule'
    Z = [m.Z for m in data_list]
    max_N = max(map(len, Z))
  elif isinstance(data_list[0], qtk.QMOutput):
    typ = 'output'
    Z = [o.molecule.Z for o in data_list if hasattr(o, 'molecule')]
    max_N = max(map(len, Z))
  else:
    qtk.warning("not supported datatype")

  if 'output' not in kwargs:
    kwargs['output'] = 'dictionary'

  xyzs = []
  Zs = []
  Es = []
  xyzStr = []

  for i in range(len(data_list)):
    if typ == 'output':
      if hasattr(data_list[i], 'molecule'):
        molecule = data_list[i].molecule
      else:
        molecule = None
      Es.append(data_list[i].Et)
    elif typ == 'molecule':
      molecule = data_list[i]

    if molecule is not None:
      zeroR = np.zeros([max_N - molecule.N, 3])
      zeroZ = np.zeros(max_N - molecule.N)
      xyz = np.vstack([molecule.R, zeroR]).tolist()
      Z = np.hstack([molecule.Z, zeroZ])
      if len(xyzs) == 0:
        #xyzs = xyz
        Zs = Z
      else:
        #xyzs = np.stack([xyzs,xyz])
        Zs = np.vstack([Zs,Z])
      xyzs.append(xyz)
  
      if kwargs['output'] == 'extended_xyz':
        xyzStr.append('%d\n' % molecule.N)
        if typ == 'output':
          xyzStr.append('%f\n' % data_list[i].Et)
        elif typ == 'molecule':
          xyzStr.append('\n')
        for I in range(molecule.N):
          r = molecule.R[I]
          xyzStr.append('%-2s % 8.4f % 8.4f % 8.4f\n'\
                        % (molecule.type_list[I], r[0], r[1], r[2]))
  xyzs = np.array(xyzs)
  
  if len(Es) > 0:
    Es = np.array(Es)
  if len(xyzStr) > 0:
    xyzStr = ''.join(xyzStr)

  if kwargs['output'] == 'dictionary':
    out = {'xyz': xyzs, 'Z': Zs}
    if len(Es) > 0:
      out['E'] = np.array(Es)
  elif kwargs['output'] == 'extended_xyz':
    out = xyzStr
  else:
    qtk.warning('not supported output format:%s' % kwargs['output'])
    out = None
  return out
