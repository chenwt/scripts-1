#!/bin/bash
#$ -cwd
#$ -l mem=4g,time=4::
#Author:J.He

bam="~/SCRATCH/projAML/WXS/callVars/PAEEYP-04A.rmdup.new.bam"
for chr in `seq 1 22`
do
  echo "MAF for ${bam} on chromosome ${chr}"
  bamName=`echo ${bam} |cut -d/ -f6 | cut -d. -f1`
  echo "PID-sampleID:"${bamName}
  echo "sh do_getMAP_aml.sh $bam $chr " | qsub -l mem=8g,time=48:: -N ${bamName}_${chr} -o ./log -e ./log -cwd
 # sh do_getMAP_aml.sh $bam $chr
done

