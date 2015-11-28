%nproc=4
%chk=/home/chang/test/
# hf gen 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine)

H

-1   1
H    0.00000000   0.00000000   0.00000000

H    0 
S   1   1.00
      0.5            1.0
****


