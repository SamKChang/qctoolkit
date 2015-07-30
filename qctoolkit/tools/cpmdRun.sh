#!/bin/bash

if [ -z "$threads" ];then
  threads=1
fi
if [ -z "$root" ];then
  root=inp
fi

inp=$1

NAME=$(echo $inp|sed 's/\.inp//g')
OUT=$(echo $inp|sed 's/inp/out/g')

mkdir $NAME
cp $root/$inp $NAME
cd $NAME

  sed -i '/MIRROR/{s/$/\n MEMORY BIG/g}' $inp

  if [ "$alchemy" = true ]; then
    maxiter=`grep -A1 $inp|tail -n 1|awk '{printf("%d",$1)}'`
    if [ "$maxiter" -gt 0 ];then
      sed -i '/MAXITER/{n;s/.*/  1/g}' $inp
    else
      sed -i '/MIRROR/{s/$/\n MAXITER\n  1/g}' $inp
    fi
    sed -i '/OPTIMIZE/{s/$/\n RESTART WAVEFUNCTION/g}' $inp
    mpirun -np $threads cpmd.x $inp > $OUT
  else
    # for cluster
    #mpirun -np $NSLOTS -mca btl tcp,self cpmd.x $inp > $OUT
    mpirun -np $threads cpmd.x $inp > $OUT
    
    if [ "$KSEg" = true ]; then
      mv RESTART.1 RESTART
      sed '/OPTIMIZE/{s/.*/ RESTART WAVEFUNCTION\n KOHN-SHAM ENERGIES\n  4\n LANCZOS PARAMETER N=5\n  50    8    20  1.D-9\n  0.05          1.D-11\n  0.01          1.D-13\n  0.0025        1.D-16\n  0.001         1.D-18/}' $inp > ks$inp
      mpirun -np $threads cpmd.x ks$inp > ks$OUT
    fi
  fi
    
  if [ "$RESTART" = true ];then
    rm GEOMETRY* KPT* LATEST*
  else
    rm RESTART* GEOMETRY* KPT* LATEST*
  fi


cd ..
