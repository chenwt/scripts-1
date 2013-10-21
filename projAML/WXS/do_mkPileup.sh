#!/bin/bash
#$ -l mem=4g,time=4::
#$ -o ./log -e ./log
#$ -cwd
#Usage: qsub do_mkPileup.sh 
#Example
#Author:J.He
#Date: Jun 20,2013

while read pid 
do 
  echo "sh mkPileup.sh ${pid}" | qsub -l mem=4g,time=8:: -o ./log -e ./log -cwd -N pileup_${pid} 
done < PID_input 
