#!/bin/bash
#$ -cwd
# Input <file with pid,each one line>
# TO DO: change the bam director and result directory each time
bamd=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/
resd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/joint1/

for chr in $(seq 1 22) X Y
do 
	if [ -f ${resd}/bams${chr}.list ]; then
		rm ${resd}/bams${chr}.list
	fi
done

while read pid; 
do 
	for chr in $(seq 1 22) X Y
	do 
		ls ${bamd}*${pid}*_split/${chr}.bam >> ${resd}/bams${chr}.list
	done
done <$1