import numpy as np
import qctoolkit as qtk
import os
import pkgutil

if pkgutil.find_loader('rdkit') is not None:
  from rdkit import Chem, rdBase
  from rdkit.Chem import AllChem, rdMolTransforms, rdMolDescriptors, Draw
  import rdkit.Chem.rdForceFieldHelpers as rcr
  from rdkit.Geometry.rdGeometry import Point3D

if pkgutil.find_loader('py3Dmol') is not None:
  import py3Dmol

def mol2rdk(mol, **kwargs):
  if 'removeHs' not in kwargs:
    kwargs['removeHs'] = False
  if type(mol) is not Chem.rdchem.Mol:
    mol.write('.tmp.xyz')
    try:
      os.system("obabel -ixyz -opdb .tmp.xyz > .tmp.pdb")
    except Exception as err:
      qtk.exit("obabel failed with error: " + str(err))
    mol = Chem.MolFromPDBFile('.tmp.pdb', **kwargs)
    os.remove('.tmp.xyz')
    os.remove('.tmp.pdb')
  return mol

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

def geopt(mol, forcefield='uff', max_iter=2000, **kwargs):

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
  return rdk2mol(rdmol)

def mol2svg(mol,
            molSize=(200,200),
            kekulize=False, 
            index=True,
            atom_label=False,
            highlight=[],
            colors={},
            sizes={},
           ):

  mol = mol2rdk(mol, removeHs=True)
    
  drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
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
  Chem.AssignAtomChiralTagsFromStructure(rdm)
  AllChem.EmbedMolecule(
    rdm,
    useExpTorsionAnglePrefs=True,
    useBasicKnowledge=True
  )
  mb = Chem.MolToMolBlock(rdm,confId=confId)
  p = py3Dmol.view(width=size[0],height=size[1])
  p.removeAllModels()
  p.addModel(mb,'sdf')
  p.setStyle({'stick':{}})
  p.setBackgroundColor(background_color)
  p.zoomTo()
  return p.show()
