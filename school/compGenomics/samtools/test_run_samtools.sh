#!/bin/bash
#Author: Jing He
#Date: 
#Last Update: 
#Example:  

##------------------------------testing code##------------------------------
#
#								this is for testing on one chromosome

##------------------------------##------------------------------
# TUMOR=$CWD'WGS/1bb5b7ef-16da-43ab-a8f4-2d98d1770e21/TCGA-AB-2978-03A-01D-0739-09_whole'
# # NORMAL=$CWD'WGS/05f440e8-6e9d-43b5-9340-5271aed310dd/TCGA-AB-2972-11A-01D-0739-09_whole.bam'
# qsub -N myjobtestX -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools view -bh tumor.bam chrX > tumor.chrX.bam"


# call somatic muation from pair of samples
qsub -N myjobCallY -e ./logs -o ./logs -b y -l mem=8G,time=2:: -S /bin/sh -cwd "samtools mpileup -DSuf wu_build36.fasta tumor.Y.bam normal.Y.bam | bcftools view -bvcgT pair -p 1.1 - > Y.var.bcf"
bcftools view Y.var.bcf | vcfutils.pl varFilter -D 60 > Y.var.vcf



##----------------filtering wiht DP4 and PV4##--

cat Y.var.vcf | grep -v "^#" | cut -f8 | sed 's/;/\t/g' | cut -f4 | sed -e 's/DP4=/\t/g'  | sed 's/\,/\t/g'| awk '{print $1+$2,'\t',$3+$4}' > ac1RefAlt.freq


# cmd="samtools mpileup -uf $REF $TUMOR.bam > $TUMOR.TEMP"
# | bcftools view -bvcg - > $TUMOR.var.raw.bcf
# # echo $cmd
# # $cmd

# # cmd="bcftools view TCGA-AB-2978-03A-01D-0739-09_whole.var.raw.bcf | vcfutils.pl varFilter -D100 > $TUMOR.flt.vcf"
# # echo $cmd
# # $cmd


# samtools view -h tumor.bam | awk '$3=="NT_113275" || /^@/' | samtools view -Sb -> onlychrNT113275.bam

