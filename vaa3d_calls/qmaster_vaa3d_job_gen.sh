#!/bin/bash

function write_vaa3d_job_config {

  outputScript=$1;
  input1=$2;
  input2=$3;
  logfile=$3;
  echo "## There are several queues available.  Please check with !ITsupport to verify which queue you should use" >> $outputScript;
  echo "#PBS -q mindscope" >> $outputScript;
  echo "# Declare that your job will use no more than some amount of memory _at_peak_" >> $outputScript;
  echo "#PBS -l vmem=16g" >> $outputScript;
  echo "# Allow up to 2 hours of walltime.  Default is 12 hours" >> $outputScript;
  echo "#PBS -l walltime=10:00:00" >> $outputScript;
  echo "# Request just one core on the host" >> $outputScript;
  echo "#PBS -l ncpus=1" >> $outputScript;
  echo "# Give your job a descriptive name. This is visible in qstat and other job reports.  Also serves as the default basename for log files" >> $outputScript;
  echo "#PBS -N ${input1}" >> $outputScript;
  echo "# Should torque automatically re-run the job on error?" >> $outputScript;
  echo "#PBS -r n" >> $outputScript;
  echo "# Merge STDOUT into STDERR" >> $outputScript;
  echo "#PBS -j oe" >> $outputScript;
  echo "# location for stderr/stdout log files _after_ job completion" >> $outputScript;
#  echo "#PBS -o ${outputScript}.out" >> $outputScript;

  echo "#" >> $outputScript;
  echo "#" >> $outputScript;

#  echo "# send email on job error" >> $outputScript;
#  echo "#PBS -m a" >> $outputScript;

  DISPLAY1=:$RANDOM;
  echo "export DISPLAY=$DISPLAY1" >> $outputScript;
  echo "Xvfb $DISPLAY1 -auth /dev/null &" >> $outputScript;
#  echo "export LD_PRELOAD=/usr/lib64/libstdc++.so.6" >> $outputScript;

  echo "export LD_LIBRARY_PATH=/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters" >> $outputScript;
  echo "/data/mat/xiaoxiaol/work/bin/bin_vaa3d_for_clusters/vaa3d -x neuron_distance -f neuron_distance -i $input1  $input2 -o $logfile " >> $outputScript;
  echo "kill %1" >> $jobScriptFile;	
}

#copy the names
input1=$1
input2=$2
logfile=$3
jobScriptFile=$4.qsub

#generate the batch script configuration
if [ -f $jobScriptFile ]; then
  rm $jobScriptFile;
fi;

write_vaa3d_job_config $jobScriptFile $input1 $input2 $logfile



