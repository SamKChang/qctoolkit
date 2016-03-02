import qctoolkit as qtk
import qctoolkit.ccs as qcs
import numpy as np
import copy, sys, re, gc, os, math
import itertools as it
import shutil

def generate(cell_x, cell_y, **kwargs):
  """
  planar graphene generator
  generate graphene using 4-atom cubic basis 
  kwargs: bond_length
  """
  if 'bond_length' in kwargs:
    bond_length = kwargs['bond_length']
  else:
    bond_length = 1.4210


  x1 = bond_length
  x2 = x1 + 0.5*bond_length
  y2 = bond_length * np.sqrt(3)/2
  x3 = x2 + bond_length
  y3 = y2
  base = np.array([
           [0.0, 0.0, 0.0],
           [ x1, 0.0, 0.0],
           [ x2,  y2, 0.0],
           [ x3,  y3, 0.0],
         ])

  graphene = qtk.Molecule()
  graphene.R = copy.deepcopy(base)

  for x in range(cell_x):
    shift_x = x * bond_length * 3
    for y in range(cell_y):
      if x>0 or y>0:
        shift_y = y * bond_length * np.sqrt(3)
        cell = base + np.array([shift_x, shift_y, 0])
        graphene.R = np.vstack(np.array([graphene.R, cell]))

  graphene.sort_coord()
  graphene.N = len(graphene.R)
  graphene.Z = [ 6 for i in range(graphene.N)]

  return graphene

def genRefInp(x, y, n_pair, **kwargs):
  """
  generate cpmd input files for graphene reference calculation
  x: number of graphene 4-atom unit cells in x
  y: number of graphene 4-atom unit cells in y
  n_pair: number of BN pairs in the system
  kwargs: max_sample
  """

  if 'path' in kwargs:
    qtk.report("creating folder", kwargs['path'])
    os.makedirs(kwargs['path'])
    os.chdir(kwargs['path'])

  print "%screating inp folder...%s"\
        % (qtk.bcolors.OKGREEN, qtk.bcolors.ENDC)
  if not os.path.exists('inp'):
    os.makedirs('inp')
  else:
    print "%sinp folder exist, baking up to back_inp...%s"\
          % (qtk.bcolors.WARNING, qtk.bcolors.ENDC)
    try:
      shutil.rmtree('back_inp')
    except:
      pass
    os.rename('inp', 'back_inp')
    os.makedirs('inp')

  N = 4*x*y

  if 'max_sample' in kwargs:
    max_sample = kwargs['max_sample']
  else:
    total = math.factorial(N-1)
    denum1 = math.factorial(n_pair)**2
    denum2 = math.factorial(N - 1 - 2*n_pair)
    max_sample = total / (denum1 * denum2)
  digit = len(str(max_sample))

  name_xyz = "gph%d-%d.xyz" % (x, y)
  name_ccs = "ccs%d-%d.txt" % (x, y)
  namev_xyz = "gph%d-%dv.xyz" % (x, y)
  namev_ccs = "ccs%d-%dv.txt" % (x, y)
  header = "gph%d%d_" % (x, y)

  ccs_file = open(name_ccs, "w")
  ccsv_file = open(namev_ccs, "w")
  print >> ccs_file, "mutation_list:\n"+\
                     " %d:%d -> 5:7\n" % (2, N) + \
                     "end"
  print >> ccsv_file, "mutation_list:\n"+\
                      " %d:%d -> 5:7\n" % (1, N-1) + \
                      "end"
  ccs_file.close()
  ccsv_file.close()

  print "%sthe following files are written: %s"\
        % (qtk.bcolors.OKGREEN, qtk.bcolors.ENDC)
  print " %s\n %s\n %s\n %s\n" % (name_xyz, namev_xyz,\
                                  name_ccs, namev_ccs)

  graphene = generate(x, y)
  graphene.write_xyz(name_xyz)
  graphenev = graphene.remove_atom(1)
  graphene.write_xyz(namev_xyz)

  space = qtk.CCS(name_xyz, name_ccs)
  space_v = qtk.CCS(namev_xyz, namev_ccs)
  flat = [item for sublist in space.mutation_list \
             for item in sublist]

  # geometry setup
  mol_base = graphene
  dx = mol_base.R[1,0]
  mol_base.sort_coord()
  dy = mol_base.R[1,1]
  center = mol_base.R[0] + [0,0,-10]
  mol_base.center(center)
  x_max = mol_base.R[-1,0]
  y_max = mol_base.R[-1,1]
  celldm = [round(x_max+dx, 4), round(y_max+dy, 4), 20, 0,0,0]


  # input setup
  gph = qtk.QMInp(name_xyz, 'cpmd', info='gph')
  gph.setCelldm(celldm)
  gph.periodic()
  gph.setSCFStep(500)
  
  gphv = qtk.QMInp(namev_xyz, 'cpmd', info='gphv')
  gphv.setCelldm(celldm)
  gphv.periodic()
  gphv.setSCFStep(500)

  name = "inp/%s%s.inp" % (header, str(0).zfill(digit))
  namev = "inp/%s%sv.inp" % (header, str(0).zfill(digit))
  gph.write(name)
  gphv.write(namev)

  c_base = [ 6 for i in range(len(flat))]

  itr = 1
  for n_bn in range(1, n_pair+1):
    for n_comb in it.combinations(range(len(flat)), n_bn):
      n_list = list(n_comb)
      rest = [index for index in range(len(flat)) \
              if index not in list(n_comb)]
      for b_comb in it.combinations(rest, n_bn):
        name = "inp/%s%s.inp" % (header, str(itr).zfill(digit))
        namev = "inp/%s%sv.inp" % (header, str(itr).zfill(digit))
        b_list = list(b_comb)
        c_list = [index for index in rest\
                  if index not in list(n_comb)]
        atom_list = copy.deepcopy(c_base)
        for i in b_list:
          atom_list[i] = 5
        for i in n_list:
          atom_list[i] = 7
        valid = True
        mutate_mol = space.generate(mutation=[atom_list])
        mutate_molv = space_v.generate(mutation=[atom_list])
        mutate_mol.find_bonds()
        if not all (key in mutate_mol.bond_types \
                    for key in ('B-B', 'N-N')):
          gph.setCenter([0,0,-celldm[2]/2])
          gph.setStructure(mutate_mol)
          gph.write(name)
          gphv.setCenter([0,0,-celldm[2]/2])
          gphv.setStructure(mutate_molv)
          gphv.write(namev)
        itr += 1
        if itr > max_sample:
          msg = "%sDONE! %d generated%s"\
                % (qtk.bcolors.OKGREEN, itr-1, qtk.bcolors.ENDC)
          sys.exit(msg)

  print "%sDONE! %d generated%s"\
         % (qtk.bcolors.OKGREEN, itr-1, qtk.bcolors.ENDC)

  itr = 1
