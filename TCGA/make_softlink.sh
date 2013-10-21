#!/bin/bash
#$ -cwd
# link the bam file downloaded to current working directory
# tab delimited file with two column required in the cwd, 
# col1:barcode col2:ananlysisid used to download data

tcgad=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/

while read line
do
	barcode=$(echo ${line}|awk '{print $1}')
	anaid=$(echo ${line} |awk '{print $2}')
	echo $barcode
	# echo $anaid
	if [ -f ${wd}${barcode}.bam ]; then 
		rm ${wd}${barcode}.bam 
	fi
	bam=$(ls ${tcgad}${anaid}/*bam |awk -F'/' '{print $NF}')
	cmd="ln -s ${tcgad}${anaid}/${bam} ${wd}${barcode}.bam"
	# echo $cmd
	$cmd
	if [ -f ${wd}${barcode}.bam.bai ]; then 
		rm ${wd}${barcode}.bam.bai
	fi
	bai=$(ls ${tcgad}${anaid}/*bai |awk -F'/' '{print $NF}')
	cmd="ln -s ${tcgad}${anaid}/${bai} ${wd}${barcode}.bam.bai"
	# echo $cmd
	$cmd
done < input_softlink.txt