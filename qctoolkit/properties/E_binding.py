import qctoolkit as qtk
#import qctoolkit.io_format.cpmd as cpmd

def E_int(mol_A, mol_B, program=qtk.setting.qmcode, **kwargs):
  mol_A = qtk.Structure(mol_A)
  mol_B = qtk.Structure(mol_B)

  inp_A = qtk.QMInp(mol_A)
  inp_B = qtk.QMInp(mol_B)

  inp_A.setVDW('DCACP')
  inp_A.write()

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
