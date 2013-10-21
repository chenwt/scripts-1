#!/bin/bash
#$ -cwd
#$ -l mem=4g,time=24::

# splitting bams into chrs

# inp2=TCGA-BH-A0B3-11B-21D-A128-09_IlluminaGA-DNASeq_whole.bam
# inp1=TCGA-BH-A0B3-01A-11D-A128-09_IlluminaGA-DNASeq_whole.bam
cnt=1
# for line in $(ls *bam)
logd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/log/
while read line 
do 
	# echo ${line}
	# echo $cnt
	cmd="~/scripts/projFocus/ceRNA/splitChr_bam.sh ${line}"
	# echo ${cmd}
	echo ${cmd} | qsub -l mem=4g,time=20:: -cwd -N split_Bam${cnt} -e ${logd} -o ${logd}  >> ${logd}qsub.log
	tail -1 ${logd}qsub.log
	let "cnt+=1"
done < input_splitChr_bam.txt