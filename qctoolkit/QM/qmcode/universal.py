import qctoolkit as qtk
import os, re

def runCode(self, parrent, name, **kwargs):
  worker, name = \
    super(parrent, self).run(kwargs['program'], name, **kwargs)

  if 'no_subfolder' not in kwargs or not kwargs['no_subfolder']:
    self.setting['root_dir'] = name

  def run():
    if 'charge' in kwargs:
      self.setChargeMultiplicity(kwargs['charge'], 1)
    inp = self.write(name, **kwargs)
    new_name = None
    if 'new_name' in kwargs:
      new_name = kwargs['new_name']
    return worker.start(inp, new_name)

  if not os.path.exists(name):
    return run()
  elif self.setting['overwrite']:
    qtk.warning("Overwrite existing folder %s" % name)
    return run()
  else:
    qtk.report("QMInp.run", "%s exists" % name)

def cornerCube(self):
  if self.setting['save_density'] or self.setting['save_wf']:
    if not self.molecule.grid:
      self.setGrid()
    if self.setting['corner_cube']:
      grid_min = [self.molecule.grid[i][0] for i in range(3)]
      self.center(grid_min)
      self.setGrid()


def PPCheck(xc, element, pp_file_str, **kwargs):
  ne = qtk.n2ve(element)
  try:
    pp_path = os.path.join(xc, element + '-q' + str(qtk.n2ve(element)))
    pp_file = os.path.join(qtk.setting.cpmd_pp_url, pp_path)
    saved_pp_path = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
    if not os.path.exists(saved_pp_path) and qtk.setting.download_pp:
      if pp_file:
        new_pp = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
        pp_content = urllib2.urlopen(pp_file).read()
        qtk.report('', 'pp file %s not found in %s. ' \
                   % (pp_file_str, qtk.setting.cpmd_pp) + \
                   'but found in cp2k page, download now...')
        new_pp_file = open(new_pp, 'w')
        new_pp_file.write(pp_content)
        new_pp_file.close()
        pp_file = new_pp
    return saved_pp_path
  except:
    qtk.warning('something wrong with pseudopotential')

def alchemyPP(xc, pp_file_str):
  pp_path = os.path.join(qtk.setting.cpmd_pp, pp_file_str)
  if not os.path.exists(pp_path):
    root, _ = os.path.splitext(pp_file_str)
    element_str = re.sub('_.*', '', root)
    fraction = float(re.sub('.*_', '', root))/100
    element1 = re.sub('2.*', '', element_str)
    element2 = re.sub('.*2', '', element_str)
    str1 = element1 + "_q" + str(qtk.n2ve(element1)) +\
           "_" + xc + '.psp'
    str2 = element2 + "_q" + str(qtk.n2ve(element2)) +\
           "_" + xc + '.psp'
    pp1 = PP(PPCheck(xc, element1, str1))
    pp2 = PP(PPCheck(xc, element2, str2))
    pp = mutatePP(pp1, pp2, fraction)
    pp.write(pp_path)

def mutatePP(pp1, pp2, fraction):
  if type(pp1) is str:
    if pp1.upper() == 'VOID':
      pp1 = PP()
    else:
      pp1 = PP(pp1)
  if type(pp2) is str:
    if pp2.upper() == 'VOID':
      pp2 = PP()
