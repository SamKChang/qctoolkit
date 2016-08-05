import numpy as np
import qctoolkit as qtk
import grid_points as gp

diag = np.diag

def getCubeGrid(inp, margin=7, step=0.2):
  coord_min = np.min(inp.mol.atom_coords(), axis=0) - margin
  coord_max = np.max(inp.mol.atom_coords(), axis=0) + margin
  steps = np.ceil((coord_max - coord_min)/float(step)) + 1
  out = np.zeros([4, 4])
  out[0, 0] = inp.molecule.N
  out[0, 1:] = coord_min
  out[1:, 0] = steps
  out[1:, 1:] = diag((coord_max - coord_min) / (steps - 1))
  return out

def getCube(inp, dv=None, cube_header=None, **kwargs):
  if dv is None:
    dv = inp.dv
  if not cube_header:
    if 'resolution' not in kwargs:
      res = 0.2
    else:
      res = float(kwargs['resolution'])
    if 'margin' not in kwargs:
      margin = 7.0
    else:
      margin = float(kwargs['resolution'])
    cube_header = getCubeGrid(inp, margin, res)

  def getSpace(i):
    coord_min = cube_header[0, 1+i]
    step = cube_header[1+i, 0]
    size = cube_header[1+i, 1+i]
    coord_max = coord_min + (step-1)*size
    return np.linspace(coord_min, coord_max, step)

  coord_axes = [getSpace(i) for i in range(3)]
  X, Y, Z = np.meshgrid(*coord_axes, indexing='ij')
  step = X.shape
  X = X.reshape(X.size)
  Y = Y.reshape(Y.size)
  Z = Z.reshape(Z.size)
  coords = np.array([X,Y,Z]).T
  rho = gp.getRho(inp, coords, new=True, dv=dv)
  rho = rho.reshape(*step)
  cube = qtk.CUBE()
  cube.build(inp.molecule, cube_header, rho)

  return cube
