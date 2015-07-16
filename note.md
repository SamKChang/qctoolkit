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
    fist switch to master branch and `git merge <branch to merge>`
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
* Default installation path: /usr/lib/python2.7/dist-packages/
<package_system>/<package_name>
**Note:** dynamic library of extensions will be compiled to 
/usr/lib/python2.7/dist-packages/<package_system>/<module_name>
* python package manager "pip" can be installed by ```sudo apt-get install python-pip```. Other packages can be installed through pip via ```sudo pip install numpy``` or upgrade via ```sudo pip install --upgrade numpy```
