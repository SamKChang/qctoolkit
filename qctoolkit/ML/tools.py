import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

from collections import OrderedDict

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

def l2_fit(X, Y, get_coefficents=False):
  """
  Y = X.dot(beta) + epsilon
  absorb constant epsilon in to X, rename beta and epsilon to w
  Y = X.dot(w)
  or
  Yi = Xi.T * wi + w0 for all i
  """
  X = np.asarray(X)
  X_base = X.copy()
  try:
    ones = np.ones(X.shape[-1])
    X = np.vstack([ones, X]).T
    w = np.linalg.solve(np.dot(X.T, X), np.dot(X.T, Y))
  except ValueError as err:
    qtk.warning('error with solve linear system: ' + str(err))
    qtk.warning('attempt transpose matrix X')
    X = X_base.T
    ones = np.ones(X.shape[-1])
    X = np.vstack([ones, X]).T
    w = np.linalg.solve(np.dot(X.T, X), np.dot(X.T, Y))
  if not get_coefficents:
    return np.dot(X, w)
  else:
    return np.dot(X, w), w


def error_measure(y, y_hat):
  y = np.asarray(y)
  err = np.abs(y-y_hat)
  err_mean = np.abs(y-y.mean())
  R = 1 - err.mean() / err_mean.mean()
  rmse = np.sqrt((err ** 2).mean())
  mae = err.mean()
  return R, rmse, mae

def coulomb_matrix(mol, n = -1, size = None, periodic_image=False,
                   sort = True, nuclear_charges = True, diagonal_power=2.4):
  if size is None:
    size = mol.N
  if size < mol.N:
    qtk.exit("matrix size too small")

  img_list = []
  if periodic_image:
    N = mol.N
    mol_ext = mol.copy()
    mol_ext.extend([3,3,3], normalize=True)
    mol = mol.copy()
    mol.R = mol_ext.R[13*N:14*N]
    for itr in range(27):
      img_list.append(mol_ext.R[itr*N:(itr+1)*N])

  positions = mol.R

  if nuclear_charges:
    charges = np.array(mol.Z)
  else:
    charges = np.ones(mol.N)
  if img_list:
    diff_list = []
    for img in img_list:
      diff = positions[:, np.newaxis, :] - img[np.newaxis, :, :]
      diff_list.append(diff)
    diff_list = np.stack(diff_list)
    diff_list = np.sqrt((diff_list ** 2).sum(axis=-1))
    
    distances = diff_list.min(axis=0)
  else:
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
  if diagonal_power or diagonal_power==0:
    charge_mask_2p4 = 0.5 * ((charges[:, np.newaxis] \
                              * charges[np.newaxis, :]) \
                              * diag_mask) ** (diagonal_power / 2.)
    cm = invR * charge_mask_Zij + charge_mask_2p4
  else:
    cm = invR * charge_mask_Zij
  if sort:
    ind = np.argsort(cm.sum(axis=-1))
    cm = cm[:, ind][ind]
  out = np.zeros([size, size])
  out[:cm.shape[0], :cm.shape[1]] = cm
  return out

def coulomb_matrices(positions, nuclear_charges = None, 
                     n = -1, sort=True):
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
             n_samples = None,
             kernels = ['laplacian'], 
             cv = None, 
             threads = 1, 
             alphas = [1e-11],
             gammas = [1e-5],
             descriptors = OrderedDict({
               coulomb_matrices: {'nuclear_charges': True}
             }),
             return_key = False,
             report = False,
            ):
  """
  return scores in the format of input parameter structure
  """

  E = data['E']

  if n_samples is None:
    n_samples = [
      int(len(E) / 10.),
      int(len(E) / 5.),
      int(len(E) / 2.),
    ]
 
  def listWrap(param):
    if '__getitem__' not in dir(param):
      param = [param]
    return param

  #descriptors = listWrap(descriptors)
  alphas = listWrap(alphas)
  gammas = listWrap(gammas)
  n_samples = listWrap(n_samples)

  #if type(descriptor_settings) is not list:
  #  descriptor_settings = [descriptor_settings]
  if not isinstance(descriptors, OrderedDict):
    if descriptors is None:
      descriptors = OrderedDict({None:None})
    elif type(descriptors) is type(coulomb_matrices):
      descriptors = OrderedDict({descriptors: {}})
    elif type(descriptors) is list \
    and type(descriptors[0]) is tuple:
      descriptors = OrderedDict(descriptors)
  if type(kernels) is not list:
    kernels = [kernels]
  
  if cv is None:
    cv = ShuffleSplit(len(E), 
                      n_iter=5,
                      test_size=.1)
  try:
    cv_fold = cv.n_iter
  except:
    cv_fold = len(cv)

  input_key = OrderedDict()
  input_key['descriptors'] = descriptors
  input_key['kernels'] = kernels
  input_key['alphas'] = alphas
  input_key['gammas'] = gammas
  input_key['n_samples'] = n_samples
  input_key['cv_fold'] = cv_fold

  output_key = OrderedDict()
  for k, v in input_key.items():
    if k == 'cv_fold':
      if cv_fold > 1:
        output_key[k] = cv_fold
    else:
      if len(v) > 1:
        output_key[k] = v

  if report:
    qtk.report("ML.tools.krrScores setting", "\n",
               "kernel:", kernels, "\n",
               "alphas:", alphas, "\n", 
               "gammas:", gammas, "\n",
               "n_samples:", n_samples, "\n",
               "cv_threads:", threads, "\n",
               "cv_fold:", cv_fold, "\n",
               "final score format: ", output_key.keys())

  all_scores = []
  for descriptor, dsetting in descriptors.items():
    descriptor_scores = []
    all_scores.append(descriptor_scores)
    if descriptor is not None:
      dsetting = copy.deepcopy(dsetting)
      if 'nuclear_charges' in dsetting\
      and dsetting['nuclear_charges']:
        dsetting['nuclear_charges'] = data['Z']
      matrix_list = descriptor(data['xyz'], **dsetting)
    else:
      matrix_list = data['X']
  
    for kernel in kernels:
      kernel_scores = []
      descriptor_scores.append(kernel_scores)
      for alpha in alphas:
        alpha_scores = []
        kernel_scores.append(alpha_scores)
        for gamma in gammas:
          gamma_scores = []
          alpha_scores.append(gamma_scores)
          kernel_ridge = KernelRidge(alpha=alpha, 
                                     gamma=gamma, 
                                     kernel=kernel)
          for n_sample in n_samples:
            if report:
              qtk.report(
                "ML.tools.krrScores, processing", "\n",
                " descriptor =",  descriptor, "\n",
                " descriptor_setting =",  dsetting, "\n",
                " kernel =",  kernel, "\n",
                " alpha =",  alpha, "\n",
                " gamma =",  gamma, "\n",
                " n_sample = ", n_sample
              )              
            cv_ = [(train[:n_sample], test) for train, test in cv]
            scores = cross_val_score(kernel_ridge, 
                                     matrix_list.reshape(
                                       len(matrix_list), -1
                                     ), 
                                     E, 
                                     cv=cv_, 
                                     n_jobs=threads, 
                                     scoring='mean_absolute_error')
            gamma_scores.append(scores)
            if report:
             qtk.report(
               "", "best score:", np.min(np.abs(scores)), "\n",
             )
  if report:             
    qtk.report("", "final format:", output_key.keys())
  if return_key:
    return np.squeeze(-np.array(all_scores)), output_key
  else:
    return np.squeeze(-np.array(all_scores))

def pack(data_array, flatten=False, **kwargs):
  list_shape = qtk.listShape(data_array)
  if type(list_shape) is int:
    return pack_list(data_array, **kwargs)
  else:
    out_shape = np.array(list_shape).shape
    out_array_dict = np.empty(out_shape, dtype=object)
    
    ind_list = []
    for ind in np.ndindex(np.array(list_shape).shape):
      ind_list.append(ind)
      data_list = np.array(data_array)[ind]
      out_array_dict[ind] = pack_list(data_list, **kwargs)
    
    if flatten:
      out_dict = {}
      keys = out_array_dict[ind].keys()
      for k in keys:
        out_dict[k] = np.concatenate(
          qtk.stack(*[out_array_dict[ind][k] for ind in ind_list])
        )
      return out_dict
    else:
      return out_array_dict

def pack_list(data_list, **kwargs):

  if isinstance(data_list[0], qtk.Molecule):
    typ = 'molecule'
    Z = [m.Z for m in data_list]
    max_N = max(map(len, Z))
  elif isinstance(data_list[0], qtk.QMOutput):
    typ = 'output'
    Z = [o.molecule.Z for o in data_list if hasattr(o, 'molecule')]
    max_N = max(map(len, Z))
  else:
    qtk.exit("not supported datatype")

  if 'output' not in kwargs:
    kwargs['output'] = 'dictionary'

  xyzs = []
  Zs = []
  Es = []
  xyzStr = []

  for i in range(len(data_list)):
    if i % 5000 == 0:
      qtk.progress("processing %d" % (i + 1))
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

def getMolecule(*args, **kwargs):
  if 'scheme' not in kwargs:
    kwargs['scheme'] = 'cheml_qmn'
  if kwargs['scheme'] == 'cheml_qmn':
    dataset = args[0]
    index = args[1]
    R_full = dataset.R[index]
    Z_full = np.atleast_2d(dataset.Z[index]).T
    Z = Z_full[Z_full > 0]
    R = R_full[(Z_full).ravel() > 0] * 0.529177249
    ZR = np.hstack([np.atleast_2d(Z).T, R])
    mol = qtk.Molecule()
    mol.build(ZR)
    return mol
