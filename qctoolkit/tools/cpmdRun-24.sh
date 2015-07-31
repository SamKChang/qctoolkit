#!/bin/bash

for inp in `ls inp`;do

  NAME=$(echo $inp|sed 's/\.inp//g')
  OUT=$(echo $inp|sed 's/inp/out/g')
  
  mkdir $NAME
  cp inp/$inp $NAME
  cd $NAME

  sed -i '/MIRROR/{s/$/\n MEMORY BIG/g}' $inp
  
  #mpirun -np $NSLOTS -mca btl tcp,self cpmd.x $inp > $OUT
  #cpmd.x $inp > $OUT &
  mpirun -np 24 cpmd.x $inp > $OUT

  if [ $2 = "KSEg" ]; then
    mv RESTART.1 RESTART
    sed '/OPTIMIZE/{s/.*/ RESTART WAVEFUNCTION\n KOHN-SHAM ENERGIES\n  4\n LANCZOS PARAMETER N=5\n  50    8    20  1.D-9\n  0.05          1.D-11\n  0.01          1.D-13\n  0.0025        1.D-16\n  0.001         1.D-18/}' $inp > ks$inp
    mpirun -np 24 cpmd.x ks$inp > ks$OUT
  fi

  rm RESTART* GEOMETRY* KPT* LATEST*
  
  cd ..
done
