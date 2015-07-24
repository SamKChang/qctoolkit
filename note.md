Git Note
========
* To sync remote repo with local one, i.e delete remote files
  `git commit -a` should be used instead of `git commit -m`
* Clean up ALL push history, set as initial commit
```bash
rm -rf .git
git init
git add .
git commit -m "Initial commit"
git remote add origin git@github.com:SamKChang/qctoolkit.git
git push -u --force origin master
```
* Branching:
  - create a branch: `git branch <new branch>`
  - switch working branch: `git checkout <target branch>`
  - merge a branch: 
    fist switch to master branch and `git merge <branch to merge>`.
    From GitHub website, Step 1: 
    From your project repository, bring in the changes and test.
`
git fetch origin
git checkout -b pandas origin/pandas
git merge master
`
    Step 2:
    Merge the changes and update on GitHub.
`
git checkout master
git merge --no-ff pandas
git push origin master
`
  - delet a branch: `git branch -d <branch to delete>`
* Remove remote files but keep the local ones:
  `git rm --cached -r <file>`
* Cloneing only branch
  `git clone -b <branch> git@address`
* Undo `git pull`, use `git reset --hard`


Python Note
===========
* Profiling: python -m cProfile -s time toAnalyse.py <args> 
from [here]
(http://stackoverflow.com/questions/582336/how-can-you-profile-a-python-script)
* Default installation path: /usr/local/lib/python2.7/dist-packages/
<package_system>/<package_name>
**Note:** dynamic library of extensions will be compiled to 
/usr/lib/python2.7/dist-packages/<package_system>/<module_name>
* python package manager "pip" can be installed by ```sudo apt-get install python-pip```. Other packages can be installed through pip via ```sudo pip install numpy``` or upgrade via ```sudo pip install --upgrade numpy```
* setup script for multipole extention can use cythonize for both C code and cython code


Python module of Open Babel
===========================
* Installation
 - compile openbabel with python\_binding enabled: ```cmake source/path/to/openbabel -DPYTHON_BINDINGS=ON```
 - install python module with pip via ```sudo pip install openbabel```

### Technical Notes
* To build molecules from atomic number and coordinates
```python
import openbabel as ob
import pybel as pb

mol = ob.OBMol() # create openbabel empty molecule object
for atom in atoms: 
  # loop through every atom's atomic number and coordinates
  new_atom = mol.NewAtom()  # empty atom object
  new_atom.SetAtomicNum(atom.Z)
  new_atom.SetVector(atom.x, atom.y, atom.z)

mol.ConnectTheDots() # find connectivity

for frag in mol.Separate():
  print pb.Molecule(frag).write("xyz")
```
**NOTE**: Memory must be properly freeed. 
Otherwise return stack/heap/segfault error
* Openbabel routine: `OBMol.DeletAtom` takes index from 1 to N
* Openbabel routine: `OBMol.DeletBond` takes index from 0 to N-1

