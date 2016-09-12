#!/bin/bash
# submission script for Sun Grid Engine cluster
# expecting a root directory, containing all jobs
# all necessary files to run a single job should 
# be wrapped by as job_directory.
# make sure exe is callable on the cluster
# syntax:
#  qsubqtk.sh exe root n_cpu 'flags' prefix

EXE=$1
ROOT=$2
NSLOTS=$3
QSLOTS=$4
FLAG=$5
if [ $# -ge 4 ];then
  PREFIX=$6
else
  PREFIX='q'
fi

EXE=`echo $EXE|sed "s/\"//g"`

if [[ $FLAG == "None" ]];then
  FLAG=''
else
  FLAG=`echo $FLAG|sed "s/'//g"`
fi

if (( $QSLOTS > 8 ));then
  paraSetup="#$ -pe orte* $QSLOTS"
else
  paraSetup="#$ -pe smp $QSLOTS"
fi

cd $ROOT
for dir in `ls -d */`; do
  cd $dir
  inp=`ls|grep -E "(inp|com|yaml)$"`
  files=`ls|grep "files$"`
  BASE=${inp%.*}
  BCHK=$BASE".chk"
  out=`echo $inp|sed 's/\.[^\.]*$/.out/g'`
  if ! [ -z "$files" ];then
    inp=$files
    out=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fi
  log=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fchk=`echo $inp|sed 's/\.[^\.]*$/.fchk/g'`
  job=$BASE$$
  cwd=$PWD

  echo "#!/bin/bash"                                   > jobsub
  echo "#$ -cwd"                                      >> jobsub
  echo "#$ -N $PREFIX${inp%.*}"                       >> jobsub
  echo "$paraSetup"                                   >> jobsub
  echo "#$ -S /bin/bash"                              >> jobsub
  if [[ $NSLOT > 8 ]];then
    echo -n "mpirun -np $NSLOTS -mca btl tcp,self "   >> jobsub
  elif [[ $NSLOT > 1 ]];then
    echo -n "mpirun -np $NSLOTS "                     >> jobsub
  fi
  echo "$EXE $inp > $out"                             >> jobsub
  if [ -z "$files" ];then
    echo "if [ -e '$log' ];then"                      >> jobsub
    echo "  mv $log $out"                             >> jobsub
    echo "fi"                                         >> jobsub
  fi
  echo "if [ -e tmp.chk ];then"                       >> jobsub
  echo "  formchk tmp.chk tmp.fchk"                   >> jobsub
  echo "  mv tmp.chk $BASE.chk"                       >> jobsub
  echo "  cubegen 1 density=scf tmp.fchk $BASE.cube"  >> jobsub
  echo "  mv tmp.fchk $BASE.fchk"                     >> jobsub
  echo "fi"                                           >> jobsub
  echo "if [ -e DENSITY ];then"                       >> jobsub
  echo "  cpmd2cube.x DENSITY"                        >> jobsub
  echo "fi"                                           >> jobsub

  sed -i "/^%nproc/{s/=.*/=$NSLOTS/g}" $inp
  sed -i "/^%chk/{s|=.*|=tmp.chk|g}" $inp
  
  qsub $FLAG jobsub
  cd ..
done
cd ..
