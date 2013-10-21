#!/bin/bash
#$ -l mem=4g,time=24
#$ -N testRe
#$ -cwd
# input: make sure all bam file in the sam directory
# output:
# example: ~/scripts/projRe/qsub_do_02_cuff.sh 

logd=/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/log/
for bam in $(ls *bam)
do 
	id=$(echo ${bam} |sed -e "s/\.bam//g")
	echo "~/scripts/projRe/do_02_cuff.sh ${bam}"|qsub -l mem=4g,time=4:: -N ${id} -e ${logd}${id} -o ${logd}${id} -cwd >> ${logd}log_qsub
	tail -1 ${logd}log_qsub
	# echo ${bam}
done
