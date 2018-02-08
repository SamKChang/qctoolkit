# Python modules for quantum chemistry applications

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

## Dependency
* [Anaconda](https://anaconda.org/)/[Miniconda](https://conda.io/miniconda.html) is recommended
* ```conda install numpy scipy cython pandas matplotlib nose pip ipython jupyter```
* ```pip install mdtraj crypto```
* ```sudo apt-get update && sudo apt-get install -y gcc g++ gfortran liblapack-dev liblapack-doc-man liblapack-doc liblapack-pic liblapack3 liblapack-test liblapacke liblapacke-dev libgsl0-dev libatlas-base-dev build-essential libffi6 libffi-dev libssl-dev libyaml-dev libpython2.7-dev python-dev freetype* libpng12-dev```

## Installation on Ubuntu 32/64 systems
* __To install__: ```cd /path/to/qctoolkit && python setup.py install --user```, ```python setup.py develop``` 
or install by pip using ```pip install qctoolkit --user```.
* Because qctoolkit depends on it self, first execution of ```import qctoolkit``` will due to the absence of the pyc files. It should work once the pyc files are created.

* __Install on Amazon Ec2__: It is tested and working on amazon Ec2 ubuntu instances. For a fresh install, all dependencies must be installed as described above.
However, it might be necessary to create temperary swap if the memory run out:
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
* __Note__ that all code are writen with __2-space indentation__.
  To change it according to pep8 standard, use the following command:
```cd /path/to/qctoolkit && find . -name "*.py"|xargs -n 1 autopep8 --in-place```
  where ```autopep8``` can be installed simply via ```pip install autopep8 --user```

## Dependent Python packages
* numpy
* scipy
* pandas
* matplotlib
* PyYAML
* cython
* psutil
* networkx
* periodictable
* mdtraj
* paramiko (newest version might be problematic, 1.17 works fine)
* And standard libraries: sys, re, os, glob, math, subprocess, multiprocessing, copy, collections, compiler.ast, shutil, fileinput, operator, inspect, xml.etree.ElementTree, ssl
* pyscf, horton, rdkit, is optional

## Implemented interfaces to QM codes
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
* Other basis:
  - [BigDFT](bigdft.org/Wiki/index.php?title=BigDFT_website)
  - [cp2k](https://www.cp2k.org/)

## Required libraries
* OpenMP
* openmpi
* gsl
(GNU Scientific Library)
* LAPACK
* libxc-3.0.0 is optional

*20150702 KYSC*
