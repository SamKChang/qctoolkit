extentions for qctoolkit
========================
Computational demanding routines can be implemented in lower level
language. Proper interface is necessary. Each code will be
compiled as run time dynamic library to be imported by python
script. The compilation is taken care of by the setup.py script.
However, proper path of dynamic library must be specified before
it can be imported. 

For each new routine, both 'setup.py' and 'c\_extention.py'
need to be updated with correct path.

