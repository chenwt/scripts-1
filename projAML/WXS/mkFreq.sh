#!/bin/bash
#$ -cwd
#$ -l mem=8g,time=10::
#$ -o ./log -e ./log
#Author:J.He

PID="PAEEYP-04A"

for i in `seq 1 22`
do
 #i=22
 awk -v chr=${i} '$1==chr' ${PID}.pileup > ${PID}_${i}.pileup
 #cmd="awk -v chr=${i} '$1==chr' ${PID}.pileup > ${PID}_${i}.pileup"
  
  #echo "Rscript getMAF_aml.R ${PID} ${i}" | qsub -l mem=8g,time=10:: -o ./log -e ./log -N Freq_${PID}_${i} -cwd 
done


#echo "Rscript getMAF_aml.R ${PID} 1" | qsub -l mem=8g,time=10:: -o ./log -e ./log -N Freq_${PID}_1 -sync yes -cwd 
#
#count=`ls *freq | wc -l `
#if [ ${count} == 22 ]; then
#
#  echo "Integrating ${count} freq files"
#  for i in `seq 1 22`
#  do
#     cat ${PID}_${i}.freq >> ${PID}.freq
#  done
#else 
#  echo ${count}" freq files generated!"
#  sleep(2m)
#fi
