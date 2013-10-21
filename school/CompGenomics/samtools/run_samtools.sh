#!/bin/bash
#Author: Jing He
#Date: Mar. 2013
#Last Update: Apr. 10, 2013
#Example: ~/scripts/school/compGenomics/samtools/run_samtools.sh

##------------------------------pipeline for call somatic mutation using samtools
# for class project Computational Genomics
# first test on pairwised sample 



TUMOR='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/WGS/79ec95e9-ab41-4fef-a673-50d2c07b3654/TCGA-AB-2966-03A-01D-0739-09_whole.bam'
NORMAL='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/WGS/82041533-b892-4039-9596-16844f2911ed/TCGA-AB-2966-11A-01D-0739-09_whole.bam'
# echo ${TUMOR:105:4} very sensitive to any folder name changing
PID=${TUMOR:105:4}
echo $PID
[ -d $PID ] || mkdir $PID
cd $PID/
echo `pwd`
# run samtools
# result folder "/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/samtools"
ln -s $TUMOR tumor.bam
ln -s $NORMAL normal.bam

# CWD='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/'
MYREF='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/samtools/wu_build36.fasta'
# # REFDBSNP='/ifs/scratch/c2b2/ys_lab/jh3283/ref/dbsnp_132_b36_resorted.vcf'
# # REFCOSMIC='/ifs/scratch/c2b2/ys_lab/jh3283/ref/b36_cosmic_v54_080711_sorted_2.vcf'

# # index the reference and bam 
#  this step only need once
# samtools faidx testref.fa 
# qsub -N myjobIR -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools faidx wu_build36.fasta"
#
mkdir logs

cmd='qsub -N myjobIT -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools index $PID/tumor.bam"'
echo $cmd  
# $cmd

cmd='qsub -N myjobIN -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools index $PID/normal.bam"'
echo $cmd
# $cmd

# splitting the bam files into different chromosomes


# for chrom in `seq 1 22` X Y
# do
# 	cmd='qsub -N myjobT$chrom -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools view -bh $TUMOR $chrom > tumor.$PID.$chrom.bam"'
# 	echo $cmd
#     # samtools view -bh tumor.bam chr${chrom} | samtools sort - chr${chrom}
#     # samtools index chr${chrom}.bam
# done

# for chrom in `seq 1 22` X Y
# do
#     cmd='qsub -N myjobT$chrom -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools view -bh $NORMAL $chrom > normal.$PID.$chrom.bam"'
# 	echo $cmd
# done


# qsub -N myjobCallY -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools mpileup -DSuf $MYREF tumor.Y.bam normal.Y.bam | bcftools view -bvcgT pair -p 1.1 - > Y.var.bcf"



