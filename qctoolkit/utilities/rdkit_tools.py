import numpy as np
import qctoolkit as qtk
import os
import pkgutil

if pkgutil.find_loader('rdkit') is not None:
  from rdkit import Chem, rdBase
  from rdkit.Chem import AllChem, rdMolTransforms, rdMolDescriptors, Draw, rdDepictor
  from rdkit.Chem.Draw import rdMolDraw2D
  import rdkit.Chem.rdForceFieldHelpers as rcr
  from rdkit.Geometry.rdGeometry import Point3D

if pkgutil.find_loader('py3Dmol') is not None:
  import py3Dmol

def smiles2mol(smiles):
  rdm = Chem.AddHs(Chem.MolFromSmiles(smiles))
  Chem.AssignAtomChiralTagsFromStructure(rdm)
  AllChem.EmbedMolecule(
    rdm,
    useExpTorsionAnglePrefs=True,
    useBasicKnowledge=True
  )
  return rdk2mol(rdm)

def mol2rdk(mol, **kwargs):
  if 'removeHs' not in kwargs:
    kwargs['removeHs'] = False
  if type(mol) is not Chem.rdchem.Mol:
    mol.write('.tmp.xyz')
    try:
      os.system("obabel -ixyz -opdb .tmp.xyz > .tmp.pdb")
    except Exception as err:
      qtk.exit("obabel failed with error: " + str(err))
    rdm = Chem.MolFromPDBFile('.tmp.pdb', **kwargs)
    os.remove('.tmp.xyz')
    os.remove('.tmp.pdb')
  return rdm

def rdk2mol(rdmol):
  conf = rdmol.GetConformer()
  ZR = []
  atoms = rdmol.GetAtoms()
  for i in range(rdmol.GetNumAtoms()):
      R = conf.GetAtomPosition(i)
      Z = atoms[i].GetAtomicNum()
      coord = [Z, R.x, R.y, R.z]
      ZR.append(coord)
  ZR = np.asarray(ZR)
  mol = qtk.Molecule()
  mol.build(ZR)
  return mol

def geopt(mol, forcefield='uff', max_iter=2000, align=False, **kwargs):

  mol = mol.copy()
  opts = mol.align(opts=True)

  if type(mol) is not Chem.rdchem.Mol:
    rdmol = mol2rdk(mol)
  else:
    rdmol = mol
  
  ff_dict = {
    'uff': AllChem.UFFOptimizeMolecule,
    'mmff94': AllChem.MMFFOptimizeMolecule
  }

  ff = ff_dict[forcefield]

  AllChem.EmbedMolecule(
    rdmol, 
    useExpTorsionAnglePrefs=True,
    useBasicKnowledge=True)
  ff(rdmol, maxIters=max_iter)
  opt_mol = rdk2mol(rdmol)
  if align:
    opt_mol.inverse_align(opts)
  return opt_mol

def forcefield_energy(mol, forcefield='uff'):
  rdm = mol2rdk(mol)
  if forcefield == 'mmff94':
      mp = AllChem.MMFFGetMoleculeProperties(rdm)
      hdl = AllChem.MMFFGetMoleculeForceField(
        rdm, mp, ignoreInterfragInteractions=False)
  elif forcefield == 'uff':
      hdl = AllChem.UFFGetMoleculeForceField(
        rdm, ignoreInterfragInteractions=False)
  return hdl.CalcEnergy()

def mol2svg(mol,
            figSize=(200,200),
            kekulize=False, 
            index=True,
            atom_label=False,
            highlight=[],
            colors={},
            sizes={},
            remove_H=False
           ):

  mol = mol2rdk(mol, removeHs=remove_H)
    
  drawer = rdMolDraw2D.MolDraw2DSVG(figSize[0],figSize[1])
  opts = drawer.drawOptions()
  kw_drawing = {}
  
  ##################
  # bare bone plot #
  ##################
  m = Chem.Mol(mol.ToBinary())
  if kekulize:
    try:
      Chem.Kekulize(m)
    except:
      m = Chem.Mol(mol.ToBinary())
  if not m.GetNumConformers():
    rdDepictor.Compute2DCoords(m)
  
  ############
  # indexing #
  ############
  if index:
    for i in range(m.GetNumAtoms()):
      opts.atomLabels[i] = m.GetAtomWithIdx(i).GetSymbol()+str(i)
          
  ##############
  # atom label #
  ##############
  if atom_label:
    for i in range(m.GetNumAtoms()):
      opts.atomLabels[i] = m.GetAtomWithIdx(i).GetSymbol()
      
  #################
  # high lighting #
  #################
  kw_hl = {
    'highlightAtoms': highlight,
    'highlightAtomColors': colors,
    'highlightAtomRadii': sizes,
    'highlightBonds': None,
  }
  kw_drawing.update(kw_hl)
  
  drawer.DrawMolecule(m, **kw_drawing)
  drawer.FinishDrawing()
  svg = drawer.GetDrawingText()
  return svg.replace('svg:','')

def show3D(mol, size=(400,400), background_color='0xeeeeee', confId=-1):
  rdm = mol2rdk(mol)
  mb = Chem.MolToMolBlock(rdm,confId=confId)
  p = py3Dmol.view(width=size[0],height=size[1])
  p.removeAllModels()
  p.addModel(mb,'sdf')
  p.setStyle({'stick':{}})
  p.setBackgroundColor(background_color)
  p.zoomTo()
  return p.show()
