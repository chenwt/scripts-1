#!/usr/bin/Rscript
#preprocessing of genotype data
#input:<snp * sample: snp genotype data>
#output: <snp * sample:after filtering>

#read into data


# 1st filtering: call_confidence < 0.1
# use a python code to do this


# 2nd filtering: HWE test. chi-square/fisher exact
# pvalue 0.001


# 3rd filtering: Minor Allele Frequency: 5% 
