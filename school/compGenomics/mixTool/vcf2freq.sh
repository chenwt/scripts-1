#!/bin/bash
#Author: Jing He
#Date: Apr. 25th
#Last Update: 
#Example:  VCF2Freq.sh
#working in each PID folder

for 
for chrom in `seq 1 22` X Y  
	do
cat | sed -e 's/##fileformat=VCFv4.1/##fileformat=VCFv4.0/g' > $chrom.somaticfiltered.vcf.temp
cat Y.var.vcf.somaticfiltered.vcf |grep -v "^#" | cut -f1,2,8 


cat Y.var.vcf.somaticfiltered.vcf | grep -v "^#"| cut -d';' -f5 | sed -e 's/DP4=//g' | sed -e 's/,/\t/g'
	done


