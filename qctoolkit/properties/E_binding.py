import qctoolkit as qtk
#import qctoolkit.io_format.cpmd as cpmd

def E_int(*segments, **kwargs):
  """
  Compute interaction energy of 1-N molecules arrangements
  a list of moleucle sould be passed as input with
  an index specifying which of the molecule is ligand
  ONLY one ligand is considered
  Note that the index is used as the final segmentation index
  """
  molecules = []
  set_ligand = False
  ligand = 0
  pocket = 1
  for inp in segments:
    if type(inp) is str:
      molecules.append(qtk.Molecule(inp))
    elif type(inp) is qtk.geometry.Molecule:
      molecules.append(inp)
    elif type(inp) is int and not set_ligand:
      set_ligand = True
      ligand = inp
      if ligand != 0: pocket = 0
    else:
      qtk.report("E_int", "wrong input format", segments)
  mol_AB = molecules[0]
  for s in range(1,len(molecules)):
    mol_AB = mol_AB + molecules[s]
  mol_AB.find_bonds()
  mol_A = mol_AB.segments[ligand]
  mol_B = molecules[pocket]
  for s in range(len(mol_AB.segments)):
    if s != ligand and s != pocket:
      mol_B = mol_B + mol_AB.segments[s]

  inp_A = qtk.QMInp(mol_A, info='interaction energy A')
  inp_B = qtk.QMInp(mol_B, info='interaction energy B')
  inp_AB = qtk.QMInp(mol_AB, info='interaction energy AB')

  inp_A.setVDW('DCACP')
  inp_A.write()
  inp_B.setVDW('DCACP')
  inp_B.write()
  inp_AB.setVDW('DCACP')
  inp_AB.write()

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
