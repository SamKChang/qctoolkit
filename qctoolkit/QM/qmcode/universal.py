def writeReturn(inp, name, **kwargs):
  if name:
    return inp, name + '.' + kwargs['extension']
  else:
    return inp

def runCode(self, parrent, name, **kwargs):
  worker, name = \
    super(parrent, self).run(kwargs['program'], name, **kwargs)

  if 'no_subfolder' not in kwargs or not kwargs['no_subfolder']:
    self.setting['root_dir'] = name
  inp, name = self.write(name, no_reset=True)
  new_name = None
  if 'new_name' in kwargs:
    new_name = self.setting['new_name']
  return worker.start(inp, new_name)

def cornerCube(self):
  if self.setting['save_density'] or self.setting['save_wf']:
    print self.molecule.R
    if not self.molecule.grid:
      self.setGrid()
    if self.setting['corner_cube']:
      grid_min = [self.molecule.grid[i][0] for i in range(3)]
      print grid_min
      self.center(grid_min)
      self.setGrid()
      print self.molecule.R
