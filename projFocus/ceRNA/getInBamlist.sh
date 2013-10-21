#!/bin/bash
#$ -cwd
# INPUT:1> PID.txt 
# 	need to change the data directory and result directory
# generate input bam list for qsub_samtools_paired_calling
# bams=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/*split/*bams
datad=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/
result=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/
# pid=$1


## for paired calling, generate list file which has the matched two file names
while read pid
do 
	if [ ! -d ${result}${pid} ]; then
		mkdir ${result}${pid}
	fi 
	cd ${result}${pid}
	if [ -f  ${result}${pid}/bams.list ]; then
		rm  ${result}${pid}/bams.list
	fi
	for line in $(ls ${datad}*${pid}*_split/*bam)
	do 
		readlink -f ${line} >> ${result}${pid}/bams.list
	done

# patients base paired wised bam.list
	for i in $(seq 1 22) X Y
	do 
		grep \/${i}.bam ${result}${pid}/bams.list| tr "\n" "\t" > ${result}${pid}/bams${i}.list
	done

done < input_getBamlist.txt

## for joint calling,