import numpy as np

class Setting(object): 
  def __init__(self): 
    """
    input setup should be manipulated here 
    to provide general interface for other program
    seperate setting is for flexibility of switching
    different input struture
    """
    # default settings
    self.theory = "PBE" 
    self.mode = "single_point" 
    self.maxstep = 1000 
    self.save_density = False 
 
    self.charge = 'auto'
    self.multiplicity = 'auto'
    self.cutoff = 100 
    self.margin = 5 
    self.center = np.array([0,0,0]) 
    self.celldm = [20,20,20,0,0,0] 
    self.lattice = np.array([[20.0,  0.0,  0.0],
                             [ 0.0, 20.0,  0.0],
                             [ 0.0,  0.0, 20.0]])
    self.unit = "Angstrom" 
    self.symmetry = "isolated" 
    self.mesh = 0 
    self.kmesh = [1,1,1] 
    self.ks_states = 0 
    self.convergence = 1.0E-5
    self.scale = [1,1,1]
    self.shift = np.array([0,0,0])
    self.vdw = 'None'

    self.set_multiplicity = False
    self.set_charge = False
    self.set_center = False
    self.set_celldm = False
    self.set_lattice = False
    self.set_margin = False
    self.set_mode = False
    self.set_step = False
    self.set_init_random = False
    self.set_scale = False
    self.set_convergence = False
    self.debug = False
    self.restart = False
    self.kpoints = False
    self.isolated = True
    self.set_shift = False
    self.set_vdw = False

  def q_symmetry(self):
    a = self.celldm[3]
    b = self.celldm[4] 
    c = self.celldm[5] 
    if self.isolated: 
      self.symmetry = 'isolated' 
      return '  ISOLATED' 
    elif a==0 and b==0 and c==0: 
      self.symmetry = 'orthorhombic' 
      return '  ORTHORHOMBIC' 
    elif a+b+c==0.5 and (a*b==0 or b*c==0 or c*a==0): 
      self.symmetry = 'triclinic' 
      return '  TRICLINIC' 

