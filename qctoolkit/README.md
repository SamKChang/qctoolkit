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

#### General I/O format

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
 - *program* # read from constructor input
 - *Ehartree* # total energy from output file in Hartree.
   By default, every unit is converted to Hartree
 - *Et* = *Ehartree* # energy to convert to other units
 - Ha2ev(): # convert to eV
 - Ha2kcal(): # convert to kcal/mol

##### analysis.py

* QMResults([path, pattern, program, parallel]) 
  # general objects contain many QM results for analysis
 - *name* # string for matplotlib data name in LaTeX math mode
 - *unit* # for energy unit conversion
 - *path* # root path of folder containing multiple QM results
   **NOTE**: back slash at the end of string MUST BE REMOVED
 - *pattern* # unix shell file pattern
 - *Et* # dictionary of {file\_name:Et}
 - *data* # numpy numerical matrix data. 
   Need to be converted from *Et*
 - sum(): # sum up all valeus of Et
 - ev(): # convert Et to eV
 - kcal(): # convert Et to kcal/mol
 - rmKey(pattern): # remove regular expression pattern from 
   key of Et
 - exKey(grouped\_pattern): # extract group(1) of grouped\_pattern
   from key of Et
 - subtract(other\_QMResults): # for every matching key
   subtract other\_QMResults.Et[key] from Et.[key]
 - subtract\_constant(constant): # linearly subtract constant
   from values of Et
 - npdata(): # convert Et to numerical numpy matrix. 
   If keys are not converted to numerical values,
   serial numbers are assigned
 - plotEt('out\_file'): # write Et to out\_file or stdout
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


### Project list:
For each project, there are specially tailored data structure.
Here is the list

#### p01\_AlGaAs:
AlGaAs prject for alchemical analysis of semiconductor crystals
