Python modules for quantum chemistry applications
=====================================================
qctoolkit is quantum chemistry tool kit. 
It meant to provide universal interface to ab initio code
to test ideas or to produce data reliably. 
The code includes Abinit, QuantumESPRESSO, Gaussian, NwChem,
CPMD, BigDFT, ... etc.
It also provide some basic molecule operations, including 
rotation, stretching, alignment, bond identification, ... etc,
and data formatting, including
xyz files, Gaussian CUBE files, V\_SIM ascii files, pdb files, ... etc.

It seems worthwile to put effort to rewrite my bash/perl/python/C 
tools in to an integrated module or package. It should boosts the
reusability, productivity, and reproducibility of my results 
generated during my PhD in Basel.
More importantly, every results should be easily reproduced, 
examined, and especially furthre developed. This package starts as 
collections of modules of format I/O, analysis, plots.
Hopefully, these modules can one day become a package for general 
purpose chemistry tool kit. 

**Installation on Ubuntu 32/64 systems**:
* __To install__: ```cd /path/to/qctoolkit && python setup.py install --user``` or install by pip using ```pip install qctoolkit --user```. 
* __Install on Amazon Ec2__: It is tested and working on amazon Ec2 ubuntu instances. For a fresh install, all dependencies must be installed
```
sudo apt-get update && sudo apt-get install -y gcc g++ gfortran liblapack-dev liblapack-doc-man liblapack-doc liblapack-pic liblapack3 liblapack-test liblapacke liblapacke-dev libgsl0-dev libatlas-base-dev build-essential libffi6 libffi-dev libssl-dev libyaml-dev libpython2.7-dev python-pip python-dev freetype* python-matplotlib cython
```
Note that matplotlib and cython are installed through ubuntu repository for convenience. It might be necessary to create temperary swap if the memory run out:
```
sudo /bin/dd if=/dev/zero of=/var/swap.1 bs=1M count=1024
sudo /sbin/mkswap /var/swap.1
sudo /sbin/swapon /var/swap.1
```
Then do pip install ```pip install qctoolkit --user```
If temerary swap is added, remove after installation:
```
sudo swapoff /var/swap.1
sudo rm /var/swap.1
```
* __To remove__:  Manually remove all created files. List of files can 
be obtained by the --record flag during install
```python setup.py install --user --record fileList.txt```
* **Note** that the ```setup.py``` script depends on python setuptools
  package. This can be installed by
```wget https://bootstrap.pypa.io/ez_setup.py --no-check-certificate -O - | python - --user```
  with superuser priviledge
* The package depends on [NumPy > 1.11.0](http://www.numpy.org/),
  [SciPy > 0.16.0](https://www.scipy.org/),
  [pandas > 0.17.1](http://pandas.pydata.org/), 
  and [matplotlib > 1.5.1](http://matplotlib.org/). 
* **Note** that newer version for many python modules are required. They must __NOT__ 
be installed via ubuntu repository. When a module is installed 
through ubuntu repository as python-modeul (e.g. python-numpy), 
import path of such module **WILL GET** highest priority. 
In other words, stable but out-dated versions will always get loaded. 
To circumvent this, 
the best solution is to use virtual enviroment and setup dependancy. 
However, it is also possible to modify the system behaviour 
by edditing the easy_install path ```/usr/local/lib/python2.7/dist-packages/easy-install.pth```
Simply comment out the second line ```/usr/lib/python2.7/dist-packages``` 
supresses the system to insert this path before PYTHONPATH.
However, for some systems this fix is still not working. 
In that case, one has to modify the python default behaviour of constructing ```sys.path```
variable. It can be done via modifying ```/usr/lib/python2.7/site.py```
and add the following two lines to the end of the file:
```
    sys.path.remove('/usr/lib/python2.7/dist-packages')
    sys.path.append('/usr/lib/python2.7/dist-packages')
```

* **Note** that all code are writen with **2-space indentation**. 
  To change it according to pep8 standard, use the following command:
```cd /path/to/qctoolkit && find . -name "*.py"|xargs -n 1 autopep8 --in-place```
  where ```autopep8``` can be installed simply via ```pip install autopep8 --user```

**Dependent Python packages**:
* numpy 1.11.0
* scipy 0.16.0
* pandas 0.17.1
* matplotlib 1.5.1
* matplotlib.pyplot
* PyYAML 3.11
* cython
* psutil
* networkx
* periodictable
* mdtraj
* paramiko (newest version might be problematic, 1.17 works fine)
* And standard libraries: sys, re, os, glob, math, subprocess, multiprocessing, copy, collections, compiler.ast, shutil, fileinput, operator, inspect, xml.etree.ElementTree
* pymol is also used for visualization
* pyscf and horton is optional

**Implemented interfaces to QM codes**:
* Gaussian basis:
  - [Gaussian](www.gaussian.com/)
  - [NWChem](www.nwchem-sw.org/index.php/Main_Page)
  - [horton](theochem.github.io/horton/)
  - [pyscf](http://sunqm.github.io/pyscf/)
* Plane wave basis:
  - [VASP](www.vasp.at)
  - [QuantumESPRESSO](www.quantum-espresso.org/)
  - [CPMD](www.cpmd.org/)
  - [Abinit](http://www.abinit.org/)
* Wavelet basis:
  - [BigDFT](bigdft.org/Wiki/index.php?title=BigDFT_website)

**Required libraries**:
* OpenMP
* openmpi
* gsl
(GNU Scientific Library)
* LAPACK
* libxc-3.0.0

*20150702 KYSC*
