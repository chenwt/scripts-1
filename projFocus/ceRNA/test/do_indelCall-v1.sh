#!/bin/bash
#By: J.He
#Desp: to qsub indelCall-v1.sh

wd='/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall'
bam='/ifs/scratch/c2b2/TCGA/data/BRCA/WXS/TCGA-A1-A0SD-01A-11D-A10Y-09.bam'
cmd="~/scripts/projFocus/ceRNA/indelCall.sh $wd $bam" 
jobName=`echo $bam|awk '{print substr($1,0,19)}'`
echo $jobName 
echo $cmd |qsub -l mem=8g,time=50:: -e $wd/log -o $wd/log -N $jobName >> $wd/log.qsub
tail -1 $wd/log.qsub



