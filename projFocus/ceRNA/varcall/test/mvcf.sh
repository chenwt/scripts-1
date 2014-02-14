#!/bin/bash
#By: J.He
#$ -cwd
#TODO: 

vcf_array=( `ls split_*bam.var.vcf.gatk.vcf ` )
VCFCONCAT=~/tools/vcftools/vcftools_current/bin/vcf-concat
echo ${vcf_array[@]}
$VCFCONCAT ${vcf_array[@]} > concatAll.vcf
