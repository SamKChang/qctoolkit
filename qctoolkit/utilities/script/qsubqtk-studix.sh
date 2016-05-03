#!/bin/bash
# submission script for Sun Grid Engine cluster
# expecting a root directory, containing all jobs
# all necessary files to run a single job should 
# be wrapped by as job_directory.
# make sure exe is callable on the cluster
# syntax:
#  qsubqtk.sh exe root n_cpu 'flags'

EXE=$1
ROOT=$2
NSLOTS=$3
FLAG=$4
PREFIX=$5

if [[ $FLAG == "None" ]];then
  FLAG=''
fi

if (( $NSLOTS > 8 ));then
  paraSetup="#$ -pe orte* $NSLOTS"
else
  paraSetup="#$ -pe smp $NSLOTS"
fi

cd $ROOT
for dir in *; do
  cd $dir
  inp=`ls|grep -E "(inp|com|yaml)$"`
  BASE=${inp%.*}
  out=`echo $inp|sed 's/\.[^\.]*$/.out/g'`
  log=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fchk=`echo $inp|sed 's/\.[^\.]*$/.fchk/g'`
  job=$BASE$$
  cwd=$PWD

  echo "#!/bin/bash"                                   > jobsub
  echo "#$ -cwd"                                      >> jobsub
  echo "#$ -N $PREFIX${inp%.*}"                       >> jobsub
  echo "$paraSetup"                                   >> jobsub
  echo "#$ -S /bin/bash"                              >> jobsub
  echo -n "mpirun -np $NSLOTS -mca btl tcp,self "     >> jobsub
  echo "$EXE $inp > $out"                             >> jobsub
  echo "if [ -e '$log' ];then"                        >> jobsub
  echo "  mv $log $out"                               >> jobsub
  echo "fi"                                           >> jobsub
  echo "if [ -e *.chk ];then"                         >> jobsub
  echo "  mv *.chk $BASE.chk"                         >> jobsub
  echo "  formchk $BASE.chk $BASE.fchk"               >> jobsub
  echo "  cubegen 1 density=scf *.fchk $BASE.cube"    >> jobsub
  echo "fi"                                           >> jobsub
  echo "if [ -e DENSITY ];then"                       >> jobsub
  echo "  cpmd2cube.x DENSITY"                        >> jobsub
  echo "fi"                                           >> jobsub

  sed -i "/^%nproc/{s/=.*/=$NSLOTS/g}" $inp
  sed -i "/^%chk/{s/=.*/=$cwd/g}" $inp
  
  qsub $FLAG jobsub
  cd ..
done
cd ..
