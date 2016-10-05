import qctoolkit as qtk
import numpy as np

def coulomb_matrix(mol, n = -2, size = 0, sort = True):
  if size == 0:
    size = mol.N
  if size < mol.N:
    qtk.exit("matrix size too small")
  positions = mol.R
  charges = mol.Z
  differences = positions[:, np.newaxis, :] \
              - positions[np.newaxis, :, :]
  distances = np.sqrt((differences ** 2).sum(axis=-1))
  distances[distances == 0] = np.nan # replace 0 for division
  invR = (distances ** n)
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

def coulomb_matrices(positions, charges, n = -2, sort=True):
  """
  return 3D numpy array of sorted Coulomb matrix
  """
  differences = positions[..., :, np.newaxis, :] \
              - positions[..., np.newaxis, :, :]
  distances = np.sqrt((differences ** 2).sum(axis=-1))
  distances[distances == 0] = np.nan
  distances[distances == 0] = -1 # replace 0 for division
  invR = (distances ** n)
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
      out['E'] = Es
  elif kwargs['output'] == 'extended_xyz':
    out = xyzStr
  else:
    qtk.warning('not supported output format:%s' % kwargs['output'])
    out = None
  return out
