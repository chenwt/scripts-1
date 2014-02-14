#!/bin/bash
#Author: Jing He
#Date: Mar. 2013
#Last Update: Feb. 10, 2014
#Example: ~/scripts/school/compGenomics/samtools/run_samtools.sh

##------------------------------pipeline for call somatic mutation using samtools
# modified  class project Computational Genomics
# first test on pairwised sample 
# CWD="/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/samtools"
sh /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh 
BASEDIR=$(dirname $(readlink -f $0))
echo $BASEDIR

FILENAME=$1
count=0
CWD=`pwd`
while read LINE
do	
	let count++
	TUMOR =`echo $LINE | awk 'BEGIN{FS="/"}{print $NF}'`      	
	$TUMOR=$CWD"/"$TUMOR
	# $NORMAL=$CWD"/"$NORMAL
	echo $TUMOR 
	# echo $NORMAL 

	# # TUMOR='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/WGS/79ec95e9-ab41-4fef-a673-50d2c07b3654/TCGA-AB-2966-03A-01D-0739-09_whole.bam'
	# # NORMAL='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/WGS/82041533-b892-4039-9596-16844f2911ed/TCGA-AB-2966-11A-01D-0739-09_whole.bam'
	# PID=${TUMOR:105:4}
	# echo $PID

	# if [[ !  -d $PID ]] ;then  mkdir $PID; fi
	# cd $BASEDIR/$PID
	# echo $(pwd)

	# if [[ !  -d logs ]] ;then  mkdir logs ; fi
	# if [ -f = "tumor.bam" ] ;then  ln -s $TUMOR tumor.bam ; fi
	# if [ -f = "normal.bam" ] ;then  ln -s $NORMAL normal.bam ; fi

	# #indexing 
	# qsub -N myjobIT -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -sync y -cwd "samtools index tumor.bam"
	# qsub -N myjobIN -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -sync y -cwd "samtools index normal.bam"

	# #splitting the bam files into different chromosomes

	# for chrom in `seq 1 22` X Y
	# do
	# 	qsub -N myjobT$chrom -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools view -bh tumor.bam $chrom > tumor.$chrom.bam"
	# done

	# for chrom in `seq 1 22` X Y
	# do
	#     qsub -N myjobN$chrom -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools view -bh normal.bam $chrom > normal.$chrom.bam"
	# done
	
	# #call raw variants
	if [ -f = "samtools_callSomaticVar.sh" ] ;then  cp ../samtools_callSomaticVar.sh . ; fi
	for chrom in `seq 1 22` X Y
	do
		qsub -N myjobCall$chrom -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "./samtools_callSomaticVar.sh $chrom"
	done

	# #call somatic variantsq
	# for chrom in `seq 1 22` X Y
	# do
	# 	ruby ~yshen/scratch/do_filter-somatic.rb -v Y.var.vcf
	# done
done < $FILENAME
echo -e "\nTotal $count pair samples calling done!"
