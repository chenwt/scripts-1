#!/bin/bash
#$ -cwd
#$ -l mem=8g,time=10::
#$ -o ./log -e ./log
#Author:J.He

PID="PAEEYP-04A"

for i in `seq 1 22`
do
  echo "mkFreq.sh ${PID} ${i}" | qsub -l mem=4g,time=10:: -o ./log -e ./log -N freq_${PID}_${i} -cwd
done

