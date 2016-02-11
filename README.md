Python 2.7 modules for quantum chemistry applications
=====================================================
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
* **Note** that the ```setup.py``` script depends on python setuptools
  package. This can be installed by
```wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python```
  with superuser priviledge
* The package depends on [NumPy](http://www.numpy.org/),
  [pandas](http://pandas.pydata.org/), 
  and [matplotlib](http://matplotlib.org/). 
* To install: ```sudo python setup.py install```
* To remove:  Manually remove all system files. List of files can 
be obtained by the --record flag during install
```sudo python setup.py install --record fileList.txt```All files

**Dependent Python packages**:
* numpy 1.9.2
* pandas 0.16.2
* matplotlib 1.4.3
* matplotlib.pyplot
* PyYAML 3.11
* cython
* psutil
* networkx
* periodictable
* mdtraj
* And standard libraries: sys, re, os, glob, math, subprocess, multiprocessing, copy, collections, compiler.ast, shutil, fileinput, operator, inspect, xml.etree.ElementTree
* pymol is also used for visualization

**Implemented interface of QM code**:
* CPMD
* VASP
* NWChem
* (horton)

**Required libraries**:
* OpenMP
* openmpi
* gsl
(GNU Scientific Library)
* LAPACK

*20150702 KYSC*
