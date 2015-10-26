#!/usr/bin/python

import qctoolkit as qtk
import qctoolkit.properties as qpty

#C = qpty.Eb('data/A.xyz', 'data/B.xyz')
C = qpty.Eb('data/AB.xyz')

#print C.write()
C.run()
