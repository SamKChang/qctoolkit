import qctoolkit as qtk
import numpy as np
import glob, os, re
from qctoolkit.MD.aimd_io import AIMDInp
from qctoolkit.MD.aimd_io import AIMDOut
from qctoolkit.MD.trajectory import xyz as xyz

class out(AIMDOut):
  def __init__(self, out_dir, **kwargs):
    AIMDOut.__init__(self, out_dir, **kwargs)
    outpath = glob.glob(os.path.join(out_dir, '*.out'))[-1]
    trjpath = os.path.join(out_dir, 'TRAJECTORY')

    outfile = open(outpath, 'r')
    out = outfile.readlines()
    outfile.close()
    pattern = re.compile('^[0-9 ]{5} {6}[A-Z][a-z ]* *[0-9-.]')
    coord_str = filter(pattern.match, out)
    cell_str = filter(lambda x: 'LATTICE VECTOR' in x, out)
    trjfile = open(trjpath, 'r')
    trjout = trjfile.readlines()
    trjfile.close()
    trj_str = filter(lambda x: 'NEW' not in x, trjout)

    data = np.array(
      [[float(r) for r in s.split( )] for s in trj_str]
    )

    self.type_list = [s.split( )[1] for s in coord_str]
    self.N = len(self.type_list)
    self.time = data[np.arange(0, len(trj_str), self.N), 0]
    self.position = data[:, 1:4].\
      reshape([len(self.time), self.N, 3])*0.529177
    self.velocity = data[:, 4:].\
      reshape([len(self.time), self.N, 3])*0.529177
    self.cell = np.array(
      [[float(c) for c in re.sub('.*:','',s).split( )]\
        for s in cell_str]
    )*0.529177 # from Bohr to Angstrom
