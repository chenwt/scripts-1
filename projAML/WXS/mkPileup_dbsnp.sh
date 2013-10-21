#!/bin/bash
#$ -cwd

#Author:J_He

#bam=$1
bam="~/SCRATCH/projAML/WXS/callVars/PAEEYP-04A.rmdup.new.bam"
bam_name=`echo ${bam} |cut -d/ -f6 | cut -d. -f1`
max_depth=10000
min_map_1=0
min_base_1=0

for i in `seq 1 22`
do
temp_bed=temp${i}.bed
cmd="samtools mpileup -DS -C 5 -d ${max_depth} -q ${min_map_q} -Q ${min_base_q} -f $REFAML -l ${temp_bed} ${bam} > ${bam_name}_${i}.pileup"
echo $cmd
echo $cmd | qsub -l mem=4g,time=4:: -o ./log -e ./log -N ${bam_name}_${i} -cwd 
done
