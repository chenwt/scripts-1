#!/bin/bash -cwd
#Author: Jing He
#Date: Apr.2nd, 2013
#Last Update: 
#Example:  


# filterVCF.sh

##------------------------------create IDs from DBSNP and COSMICs
cat /ifs/scratch/c2b2/ys_lab/jh3283/ref/dbsnp_132.b36.vcf | grep -v "^#" | cut -f3 > rsID_dbsnp_132_b36.txt
cat /ifs/scratch/c2b2/ys_lab/jh3283/ref/b36_cosmic_v54_080711_sorted_2.vcf | grep -v "^#" | cut -f3 > rsID_cosmic_v54.txt

--exclude <filename>
