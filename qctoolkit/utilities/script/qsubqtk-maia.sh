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
if [ $# -ge 5 ];then
  PREFIX=$6
else
  PREFIX='q'
fi

CWD=`pwd`

EXE=`echo $EXE|sed "s/\"//g"`

if [[ $FLAG == "None" ]];then
  FLAG=''
else
  FLAG=`echo $FLAG|sed "s/'//g"`
fi

if (( $QSLOTS > 8 ));then
  paraSetup="#$ -pe mpi $QSLOTS"
else
  paraSetup="#$ -pe smp $QSLOTS"
fi

cd $ROOT
# submit job in each directory under ROOT
for dir in `ls -d */`; do
  cd $dir
  mydir=$PWD
  inp=`ls|grep -E "(inp|com|yaml)$"`
  BASE=${inp%.*}
  BCHK=$BASE".chk"
  out=`echo $inp|sed 's/\.[^\.]*$/.out/g'`
  files=`ls|grep "files$"`
  if ! [ -z "$files" ];then
    sed -i '1iiomode 1' $inp
    inp=$files
    out=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fi
  log=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fchk=`echo $inp|sed 's/\.[^\.]*$/.fchk/g'`
  job=$BASE$$
  cwd=$PWD

  # job setup
  echo "#!/bin/bash"                                   > jobsub
  echo "module load OpenMPI/1.6.4-GCC-4.7.2"          >> jobsub
  echo "module load hwloc/1.6.2-GCC-4.7.2"            >> jobsub
  echo "module load ABINIT/8.0.8-goolf-1.7.20"        >> jobsub
  echo "module load gaussian"                         >> jobsub
  echo "#$ -cwd"                                      >> jobsub
  echo "#$ -N $PREFIX${inp%.*}"                       >> jobsub
#  echo 'NCPU=`nproc --all`'                           >> jobsub
#  echo -n "if [[ $QSLOTS > "                          >> jobsub
#  echo '$NCPU ]];then'                                >> jobsub
#  echo "  #$ -pe mpi $QSLOTS"                         >> jobsub
#  echo "else"                                         >> jobsub
#  echo "  #$ -pe smp $QSLOTS"                         >> jobsub
#  echo "fi"                                           >> jobsub
  echo "$paraSetup"                                   >> jobsub
  echo "#$ -S /bin/bash"                              >> jobsub

  # construct scratch folder
  echo "rm -rf /tmp/$USER/$job"                       >> jobsub
  echo "mkdir -p /tmp/$USER"                          >> jobsub
  echo "mkdir /tmp/$USER/$job"                        >> jobsub
  echo "cp -r * /tmp/$USER/$job"                      >> jobsub
  echo "cd /tmp/$USER/$job"                           >> jobsub

  # mpi jobs
  echo 'NCPU=`nproc --all`'                           >> jobsub
  echo -n "if [[ $NSLOTS > "                          >> jobsub
  echo '$NCPU ]];then'                                >> jobsub
  echo -n "  mpirun -np $NSLOTS -mca btl tcp,self "   >> jobsub
  echo "$EXE $inp > $out"                             >> jobsub
  echo "elif [[ $NSLOTS > 1 ]];then"                  >> jobsub
  echo "  mpirun -np $NSLOTS $EXE $inp > $out"        >> jobsub
  echo "else"                                         >> jobsub
  echo "  $EXE $inp > $out"                           >> jobsub
  echo "fi"                                           >> jobsub

  # unified outputs
  # gaussian output
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
  # cpmd density file
  echo "if [ -e DENSITY ];then"                       >> jobsub
  echo "  cpmd2cube.x DENSITY"                        >> jobsub
  echo "fi"                                           >> jobsub

  # cleanup scratch
  echo "cd .."                                        >> jobsub
  echo "cp -rT /tmp/$USER/$job/* $mydir"              >> jobsub
  echo "cd $mydir"                                    >> jobsub
  echo "rm -rf /tmp/$USER/$job"                       >> jobsub

  sed -i "/^%nproc/{s/=.*/=$NSLOTS/g}" $inp
  sed -i "/^%chk/{s|=.*|=$PWD/$BCHK|g}" $inp
  
  if [[ $EXE == "g09" ]];then
    if [[ $FLAG == *"runtime"* ]];then
      walltime=`echo $FLAG|sed 's|.*runtime *=\([0-9:]*\).*|\1|g'`
    fi
    if [ ${walltime+x} ];then
      cmd="gsub $inp --runtime=$walltime"
    else
      cmd="gsub $inp"
    fi
    eval $cmd
    echo "# $cmd" >> *.sh
    echo "formchk $PWD/$BCHK" >> *.sh
  else
    qsub $FLAG jobsub
  fi
  cd ..
done
cd ..
