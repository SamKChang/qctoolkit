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
  inp=`ls|grep $dir`
  BASE=${f%.*}
  out=`echo $inp|sed 's/\.[^\.]*$/.out/g'`
  job=$BASE$$
  cwd=$PWD

  echo "#!/bin/bash"                                   > jobsub
  echo "#$ -cwd"                                      >> jobsub
  echo "#$ -N q${inp%.*}"                             >> jobsub
  echo "$paraSetup"                                   >> jobsub
  echo "#$ -S /bin/bash"                              >> jobsub
  echo -n "mpirun -np $NSLOTS -mca btl tcp,self "     >> jobsub
  echo "$EXE $inp > $out"                             >> jobsub
  
  qsub $FLAG jobsub
  cd ..
done
