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
FLAG=$4
if [ $# -ge 4 ];then
  PREFIX=$5
else
  PREFIX='q'
fi

if [[ $FLAG == "None" ]];then
  FLAG=''
else
  FLAG=`echo $FLAG|sed "s/'//g"`
fi

if (( $NSLOTS > 8 ));then
  paraSetup="#$ -pe orte* $NSLOTS"
else
  paraSetup="#$ -pe smp $NSLOTS"
fi

cd $ROOT
for dir in *; do
  cd $dir
  mydir=$PWD
  inp=`ls|grep -E "(inp|com|yaml)$"`
  BASE=${inp%.*}
  out=`echo $inp|sed 's/\.[^\.]*$/.out/g'`
  log=`echo $inp|sed 's/\.[^\.]*$/.log/g'`
  fchk=`echo $inp|sed 's/\.[^\.]*$/.fchk/g'`
  job=$BASE$$
  cwd=$PWD

  echo "#!/bin/bash"                                   > jobsub
  echo "module load openmpi/gnu/1.6.5"                >> jobsub
  echo "module load gaussian"                         >> jobsub
  echo "#$ -cwd"                                      >> jobsub
  echo "#$ -N $PREFIX${inp%.*}"                       >> jobsub
  echo "$paraSetup"                                   >> jobsub
  echo "#$ -S /bin/bash"                              >> jobsub
  echo "rm -rf /tmp/$USER/$job"                       >> jobsub
  echo "mkdir -p /tmp/$USER"                          >> jobsub
  echo "mkdir /tmp/$USER/$job"                        >> jobsub
  echo "cp -r * /tmp/$USER/$job"                      >> jobsub
  echo "cd /tmp/$USER/$job"                           >> jobsub
  echo -n "mpirun -np $NSLOTS -mca btl tcp,self "     >> jobsub
  echo "$EXE $inp > $out"                             >> jobsub
  echo "if [ -e '$log' ];then"                        >> jobsub
  echo "  mv $log $out"                               >> jobsub
  echo "fi"                                           >> jobsub
  echo "for chk in *.chk;do"                          >> jobsub
  echo "  if [ -e $f ];then"                          >> jobsub
  echo "    cp $f $BASE.chk"                          >> jobsub
  echo "  fi"                                         >> jobsub
  echo "done"                                         >> jobsub
  echo "if [ -e $BASE.fchk ];then"                    >> jobsub
  echo "  formchk $BASE.chk $BASE.fchk"               >> jobsub
  echo "  cubegen 1 density=scf *.fchk $BASE.cube"    >> jobsub
  echo "fi"                                           >> jobsub
  echo "if [ -e DENSITY ];then"                       >> jobsub
  echo "  cpmd2cube.x DENSITY"                        >> jobsub
  echo "fi"                                           >> jobsub
  echo "cd .."                                        >> jobsub
  echo "cp -rT /tmp/$USER/$job $mydir"                >> jobsub
  echo "cd $mydir"                                    >> jobsub
  echo "rm -rf /tmp/$USER/$job"                       >> jobsub

  sed -i "/^%nproc/{s/=.*/=$NSLOTS/g}" $inp
  sed -i "/^%chk/{s|=.*/=$cwd|g}" $inp
  
  qsub $FLAG jobsub
  cd ..
done
cd ..
