import qctoolkit as qtk
import re
#import qctoolkit.io_format.cpmd as cpmd

def EbRun(EbObject, **kwargs):
  """
  general rapper for parallel execution
  """
  return EbObject.run(**kwargs)

class Eb(object):
  """
  Compute interaction energy of 1-N molecules arrangements
  a list of moleucle sould be passed as input with
  an index specifying which of the molecule is ligand
  ONLY one ligand is considered
  Note that the index is used as the final segmentation index
  """
  def __init__(self, *segments, **kwargs):
    molecules = []
    set_ligand = False
  
    # flexible input format
    ligand = 0
    pocket = 1
    for inp in segments:
      if type(inp) is str:
        molecules.append(qtk.Molecule(inp))
      elif type(inp) is qtk.geometry.Molecule:
        molecules.append(inp)
      elif type(inp) is int and not set_ligand:
        if inp <= 0:
          qtk.exit("ligand starts from 1")
        set_ligand = True
        qtk.report("property.E_binding", "set ligand", inp - 1)
        ligand = inp - 1
        if ligand != 0: pocket = 0
      else:
        qtk.report("Eb", "wrong input format", segments)

    if len(segments) == 1 and type(segments[0]) is str:
      self.header = qtk.fileStrip(segments[0]) + "-"
    else:
      self.header = ''
  
    # construct total system AB
    self.mol_AB = molecules[0]
    for s in range(1,len(molecules)):
      self.mol_AB = self.mol_AB + molecules[s]
    # construct ligand A and pocket B
    if 'margin' not in kwargs:
      kwargs['margin'] = 5
    self.celldm = self.mol_AB.setMargin(kwargs['margin'])
    self.mol_AB.find_bonds()
    self.mol_A = self.mol_AB.segments[ligand]
    self.mol_B = self.mol_AB.segments[pocket]
    for s in range(len(self.mol_AB.segments)):
      if s != ligand and s != pocket:
        self.mol_B = self.mol_B + self.mol_AB.segments[s]
  
    if 'program' not in kwargs:
      kwargs['program'] = qtk.setting.qmcode
    if 'vdw' not in kwargs:
      if kwargs['program'] == 'cpmd':
        kwargs['vdw'] = 'DCACP'
      elif kwargs['program'] == 'vasp':
        kwargs['vdw'] = 'mbd_iter'

    self.kwargs = kwargs
    self.setInp()  

  def setInp(self):
    kwargs = self.kwargs
    print kwargs
    kwargs['celldm'] = self.celldm
    kwargs['info'] = self.header + 'Eb_A'
    self.A = qtk.QMInp(self.mol_A, **kwargs)
    kwargs['info'] = self.header + 'Eb_B'
    self.B = qtk.QMInp(self.mol_B, **kwargs)
    kwargs['info'] = self.header + 'Eb_AB'
    self.AB = qtk.QMInp(self.mol_AB, **kwargs)

  def view(self):
    self.A.view('A')
    self.B.view('B')
    self.AB.view('AB')

  def run(self, **kwargs):
    header = re.sub('-', '', self.header)
    E_A = self.A.run('Eb_A-' + header, **kwargs)
    E_B = self.B.run('Eb_B-' + header, **kwargs)
    E_AB = self.AB.run('Eb_AB-' + header, **kwargs)

    Eb = E_AB - E_A - E_B
    return Eb





#  inp_A.write()
#  inp_B.write()
#  inp_AB.write()














#  if program == 'cpmd':
#
#    if 'vdw' in kwargs:
#      vdw = kwargs['vdw']
#    else:
#      vdw = 'DCACP'
#
#    inp_A = cpmd.inp(mol_A, 
#                     '2 body interaction energy part A', 
#                     **kwargs)
#    inp_B = cpmd.inp(mol_B, 
#                     '2 body interaction energy part B', 
#                     **kwargs)
#
#    if vdw == 'DCACP':
#      inp_A.setting
#
#    print inp_A.structure.type_list
