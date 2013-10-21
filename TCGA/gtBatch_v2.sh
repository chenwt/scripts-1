#!/bin/bash 
#J.HE
#Input: <tab delimted file, barcode | analysis id>
#Ouput: downloaded bam file, and splitted bam file
#PURPOSE: 1. downlaod bam 
		 # 2. move bam to current folder and rename bam 
		  # 3, trigger splitting bam by chr
#

logd='log/'
FILENAME=$1
count=0
#yslab key
#KEY="/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/key/cghub.key"
#aclab key
KEY="/ifs/scratch/c2b2/TCGA/soft/GeneTorrent/mykey.pem"
while read LINE
do	
	anaid=$(echo ${LINE} | awk '{print $2}')
	barcode=$(echo ${LINE} | awk '{print $1}')
	let count++
      	echo "No. $count file anaid ${anaid} downloaded" >> ${logd}logs.downloaded
      	cmd="gtdownload -v -c $KEY -d ${anaid}"
	# echo $cmd 
	$cmd >> logs.download 
	echo "$LINE" >> logs.downloaded

#  renmae file
	# cmd="mv ${anaid}.bam ${barcode}.bam"
	cmd="mv ${anaid}/*bam ${barcode}.bam"
	# echo $cmd 
	$cmd
	cmd="mv ${anaid}/*bam.bai ${barcode}.bam.bai"
	# cmd="mv ${anaid}.bam.bai ${barcode}.bam.bai"
	# echo $cmd 
	$cmd

	cmd="~/scripts/projFocus/ceRNA/splitChr_bam.sh ${barcode}.bam"
	# scd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/"
	# echo ${cmd}
	echo ${cmd} | qsub -l mem=4g,time=20:: -cwd -N split_${barcode} -e ${logd} -o ${logd}  >> ${logd}qsub.log
	tail -1 ${logd}qsub.log

done < $FILENAME
echo -e "#----------------END------------------"
