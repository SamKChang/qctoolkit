import qctoolkit as qtk
import numpy as np

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
