%nproc=4
%chk=/home/chang/test/
# hf gen 6d 10f nosymm Scf(maxcycle=1000,verytight) int(grid=ultrafine)

H

-1   1
H    0.00000000   0.00000000   0.00000000

H    0 
S   1   1.00
      0.5            1.0
S   3   1.00
     82.6400000              0.0020060        
     12.4100000              0.0153430        
      2.8240000              0.0755790
****


