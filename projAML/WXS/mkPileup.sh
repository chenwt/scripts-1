#!/bin/bash
#$ -cwd
#$ -l mem=4g,time=4:: -o ./log -e ./log 
#$ -N ${bam_name}
#Author:J_He
#Usage: $0 <PID-SAMPLECODE>; exmaple: $0 PAEEYP-04A

date=$(date)
echo $date >> $(basename $0).log

bam="~/SCRATCH/projAML/WXS/callVars/"$1".rmdup.new.bam"
bam_name=`echo ${bam} |cut -d/ -f6 | cut -d. -f1`
max_depth=10000
min_map_1=0
min_base_1=0



temp_bed="/ifs/scratch/c2b2/ac_lab/jh3283/ref/Exome_Targeted_Regions.BED"
cmd="samtools mpileup -DS -C 5 -d ${max_depth} -q ${min_map_q} -Q ${min_base_q} -f $REFAML -l ${temp_bed} ${bam} > ${bam_name}.pileup"
echo $cmd
echo $cmd | qsub -l mem=4g,time=10:: -o ./log -e ./log -N ${bam_name} -cwd>> $(basename $0).log

cat $(basename $0).log
