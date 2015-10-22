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
        molecules.append(qtk.Molecule(inp, **kwargs))
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
      self.header = qtk.fileStrip(segments[0])
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
    self.mol_AB.find_bonds(**kwargs)
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
    kwargs['celldm'] = self.celldm
    kwargs['info'] = 'Eb_A-' + self.header
    self.A = qtk.QMInp(self.mol_A, **kwargs)
    kwargs['info'] = 'Eb_B-' + self.header
    self.B = qtk.QMInp(self.mol_B, **kwargs)
    kwargs['info'] = 'Eb_AB-' + self.header
    self.AB = qtk.QMInp(self.mol_AB, **kwargs)

  def view(self):
    self.A.view('A')
    self.B.view('B')
    self.AB.view('AB')

  def run(self, **kwargs):
    E_A = self.A.run('Eb_A-' + self.header, **kwargs)
    E_B = self.B.run('Eb_B-' + self.header, **kwargs)
    E_AB = self.AB.run('Eb_AB-' + self.header, **kwargs)

    Eb = E_AB - E_A - E_B
    return Eb

  def write(self, *args, **kwargs):
    if len(*args) > 0:
      self.A.write('Eb_A-'+args[0], **kwargs)
      self.B.write('Eb_B-'+args[0], **kwargs)
      self.AB.write('Eb_AB-'+args[0], **kwargs)
    else:
      self.A.write('Eb_A'+self.header, **kwargs)
      self.B.write('Eb_B-'+self.header, **kwargs)
      self.AB.write('Eb_AB-'+self.header, **kwargs)
    



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
