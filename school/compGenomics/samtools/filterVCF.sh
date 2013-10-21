#!/bin/bash -cwd
#Author: Jing He
#Date: Apr.2nd, 2013
#Last Update:Apr.14,2013 
#Example: FilterVCF.sh Y.var.vcf 


# filterVCF.sh

##------------------------------create IDs from DBSNP and COSMICs
#cat /ifs/scratch/c2b2/ys_lab/jh3283/ref/dbsnp_132.b36.vcf | grep -v "^#" | cut -f3 > rsID_dbsnp_132_b36.txt
#cat /ifs/scratch/c2b2/ys_lab/jh3283/ref/b36_cosmic_v54_080711_sorted_2.vcf | grep -v "^#" | cut -f3 > rsID_cosmic_v54.txt

#--exclude <filename>

#extracting DP4 and PV4 
cat Y.var.vcf | grep -v "^#" | cut -f1,2,8 | sed '/INDEL/d' | awk '/DP4/ && /PV4/{print}' | awk -F'[\t;]' '{print $1,$2,$6,$7,$9,$10}' >temp.1.vcf

cat temp.1.vcf | sed -e 's/AC1=\([0-9]*\)//g' | sed -e 's/MQ=\([0-9]*\)//g' | sed -e 's/FQ=\([-]*[0-9]*\.*[0-9]*\)//g' | awk '{print $1"\t"$2"\t"$3"\t"$4}' > temp.2.vcf 

#filetering with DP4, altallel > 4, alt_Allel_rev=alt_allel_fwd > 2, PV4 >

cat temp.2.vcf | sed -e 's/\(DP4=\)//g' |sed -e 's/\(PV4=\)//g' | awk -F'[, ]' '{if($6>2&&$7>2) print $0}' > temp.4.vcf


