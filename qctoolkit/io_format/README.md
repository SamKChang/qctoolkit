I/O format of QM packages
-------------------------

This folder contains I/O formats of some commonly used QM codes.
New format can be implemented.

##### setting.py

* QMSetting() # general settings
 - *theory* = "PBE"
 - *mode* = "single\_point"
 - *maxstep* = 10000
 - *save\_density* = False
* QMSetting() # plane wave settings
 - *cutoff* = 100 
 - *margin* = 5 
 - *center* = np.array([0, 0, 0]) 
 - *celldm* = [20, 20, 20, 0, 0, 0]
 - *unit* = "Angstrom"
 - *symmetry* = "isolated"
 - *mesh* = 0 
 - *kmesh* = [1,1,1]
 - *ks\_states* = 0

##### cpmd.py

* inp(structureFile, info) # info field is requied
 - *setting* = QMSetting()
 - *set\_center* = False
 - *set\_celldm* = False
 - *set\_margin* = False
 - *set\_mode* = False
 - *debug* = False
 - *restart* = False
 - PPStringDefault(atom\_type):
 - write(file\_name):
* out(cpmd\_outFile)
