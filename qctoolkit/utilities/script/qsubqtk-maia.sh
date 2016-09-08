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
# submit job in each directory under ROOT
for dir in `ls -d */`; do
  cd $dir
  mydir=$PWD
  inp=`ls|grep -E "(inp|com|yaml)$"`
  BASE=${inp%.*}
  BCHK=$BASE".chk"
  out=`echo $inp|sed 's/\.[^\.]*$/.out/g'`
  log=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fchk=`echo $inp|sed 's/\.[^\.]*$/.fchk/g'`
  job=$BASE$$
  cwd=$PWD

  # job setup
  echo "#!/bin/bash"                                   > jobsub
  echo "module load OpenMPI/1.8.4-GCC-4.8.4"          >> jobsub
  echo "module load gaussian"                         >> jobsub
  echo "#$ -cwd"                                      >> jobsub
  echo "#$ -N $PREFIX${inp%.*}"                       >> jobsub
  echo "$paraSetup"                                   >> jobsub
  echo "#$ -S /bin/bash"                              >> jobsub

  # construct scratch folder
  echo "rm -rf /tmp/$USER/$job"                       >> jobsub
  echo "mkdir -p /tmp/$USER"                          >> jobsub
  echo "mkdir /tmp/$USER/$job"                        >> jobsub
  echo "cp -r * /tmp/$USER/$job"                      >> jobsub
  echo "cd /tmp/$USER/$job"                           >> jobsub

  # mpi jobs
  echo -n "mpirun -np $NSLOTS -mca btl tcp,self "     >> jobsub
  echo "$EXE $inp > $out"                             >> jobsub

  # unified outputs
  # gaussian output
  echo "if [ -e '$log' ];then"                        >> jobsub
  echo "  mv $log $out"                               >> jobsub
  echo "fi"                                           >> jobsub
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
  echo "cp -rT /tmp/$USER/$job $mydir"                >> jobsub
  echo "cd $mydir"                                    >> jobsub
  echo "rm -rf /tmp/$USER/$job"                       >> jobsub

  sed -i "/^%nproc/{s/=.*/=$NSLOTS/g}" $inp
  sed -i "/^%chk/{s|=.*|=tmp.chk|g}" $inp
  
  qsub $FLAG jobsub
  cd ..
done
cd ..
