#!/usr/bin/python

import os
import sys
import re

if os.path.exists(sys.argv[1]):

  module_list = []
  module = re.compile("[a-zA-Z0-9_]*\.py")
  for root, dirs, files in os.walk(sys.argv[1]):
    for x in files:
      if (not re.match("__init__", x)):
        print x
        if re.match(module, x):
          module_list.append(re.sub("\.py", "", x))
        else:
          print "not matched!"
  print module_list

  os.chdir(sys.argv[1])
  init = open("__init__.py", "w")
  print >> init, "__all__ =", module_list
  init.close()
  
    #print "root: ", root
    #print "dirs: ", dirs
    #print "files: ", files
  
else:
  print "project: '", sys.argv[1], "' not exist!"
