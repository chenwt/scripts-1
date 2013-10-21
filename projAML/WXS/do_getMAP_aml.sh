#!/bin/bash
#$ -cwd
#$ -o ./log -j yes
#$ -l mem=4g,time=4::

#Author:J.He
#Date: Jun 14,2013
#Last Update: 

bam=$1
#bam="/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/callVars/PAEEYP-04A.rmdup.new.bam"
chr=$2
Rscript getMAF_aml.R $bam $chr

