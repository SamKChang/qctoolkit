import numpy as np

# parameters for plane wave QM calculations 
# includeing CPMD, espresso, NwChem 
class QMSetting(object): 
  def __init__(self): 
 
    # general default settings 
    self.theory = "PBE" 
    self.mode = "single_point" 
    self.maxstep = 1000 
    self.save_density = False 
 
    # plane wave default settings 
    self.cutoff = 100 
    self.margin = 5 
    self.center = np.array([0,0,0]) 
    self.celldm = [20,20,20,0,0,0] 
    self.unit = "Angstrom" 
    self.symmetry = "isolated" 
    self.mesh = 0 
    self.kmesh = [1,1,1] 
    self.ks_states = 0 

