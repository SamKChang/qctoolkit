Is it editable via DroidEdit?

Manual
======
**CAUTION**: No foolproof check is implemented! 
Proper file formats, data structures, I/O interfaces are assumed. 

### Ideal work flow:
Ideally, every input files will be enumerated from very simple
"template" file using a simple `geninp.py` script. 
With predefined file/data structure, every input
file should be in "inp" folder under some "root" folder 
with proper submission script for cluster. 
For each input file in "inp" folder, the submission script
should generate independent folder for calculation on the cluster.
The results contained in "root" folder should be 
analyzed and plotted by `genPlot.py` script.
```
genInput.py <template> <root>
qsub submit.sh
genPlot.py <root>
```
Note that special file structrue can be added in project folder.

### Module/Class list:
This section is a brief reference to list of modules, classes, 
attributes, and methods. **Module** are highlighted by bold
font with expected argument for constructor, while
classes are enumerated by lists. *Attributes* of each
class are denoted by italic with default values. 
Methods are followed by semicolon.
Short description of the class/attribute/method are followed
by a hash symbol as a convention for comments.


### qctoolkit - General I/O format

##### geometry.py

* Molecule()
 - *N* = 0
 - *R* = np.atleast\_2d(np.array([]))
 - *type\_list* = 'None'
 - *Z* = 0
 - *charge* = 0
 - *index* = 0
 - sort(): # sort atom types for counting, used for plane wave code
 - center(center\_coord): 
   # shift molecules to center\_coord as (0, 0, 0)
 - read\_xyz('xyz\_file'): # read from xyz\_file
 - print\_xyz(): # Print to stdout
 - write\_xyz('out\_file'): # write to out\_file
 - read\_cpmdinp('cpmd\_inp'): # read from CPMD input file

##### qmio.py

* QMInp(structur\_inp, program, info) 
  # general interface for various QM code input
 - *program* # read from constructor input. 
   Used to determine input format
 - *info* # read from constructor input vector.
   Used for result analysis.
 - setAtom(atom\_list, atom\_string): 
   # manually set atom symbols for PP or charges
 - setCorner(corner\_coord): # set corner of periodic box 
   for plan wave code
 - setCenter(center\_coord): # linearly move structure to 
   selected center
 - setCelldm(celldm): # periodic box for plan wave code
 - setMargin(margin): # distance to periodic box
 - setMode(mode): # e.g. single\_point, geopt
 - setChargeMultiplicity(charge, multiplicity):
 - setTheory(theory): # e.g pbe, pbe0, blyp, b3lyp
 - setSCFStep(step):
 - restart():
 - debug():
* QMOut(output\_file, program)
  # general interface for various QM code output
 - *program* # read from constructor input
 - *Ehartree* # total energy from output file in Hartree.
   By default, every unit is converted to Hartree
 - *Et* = *Ehartree* # energy to convert to other units
 - Ha2ev(): # convert to eV
 - Ha2kcal(): # convert to kcal/mol

##### analysis.py

* QMResults([path, pattern, program, parallel]) 
  # general objects contain many QM results for analysis
 - *path* # root path of folder containing multiple QM results
   **NOTE**: back slash at the end of string MUST BE REMOVED
 - *pattern* # unix shell file pattern
* QMData([path, pattern, program, parallel])
  # pandas DataFrame wrapper
 - ev(): # convert Et to eV
 - kcal(): # convert Et to kcal/mol
* ScatterPlot(QMResults\_predicted, QMResults\_true)
 - *pred* # QMResults\_predicted
 - *true* # QMResults\_true
 - *data* # numpy matrix for plotting
 - plot(): # generate default matplotlib object
 - show(): # show default plot on stdout
 - save(out.pdf): save default plot to out.pdf

##### utilities.py

 - R(theta, u): # rotation matrix
 - n2ve('atom\_type'): # valence electron
 - Z2n(Z): # atomic number to atom\_type
 - n2Z('atom\_type'): # atom\_type to atomic number

### qctoolkit/io\_format - QM interface

##### setting.py

* QMSetting

##### cpmd.py

* inp
* out

### qctoolkit/src - Pure C extension

### Project list:
For each project, there are specially tailored data structure.
Here is the list

#### qctoolkit/projects/Basel/p01\_AlGaAs:
AlGaAs prject for alchemical analysis of semiconductor crystals
