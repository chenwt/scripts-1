#!/bin/bash
#$ -cwd
# get_Pvaluefor_consistent.sh
# TODO: 


##----------------------------
#step 1 calculate mutations number
wd1=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result/
# wd2=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/muTect/varAccepted/
wd2=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/muTect/result/
pidfile=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/PID_16.txt
echo 'pid samtools_tn samtools_rn gatk_tn\gatk_rn\n' > temp_mutation_count
while read pid 
do
	# filtered
	# vcf1=${pid}-TuA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv.vcf
	# vcf2=${pid}-ReA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv.vcf
	# raw
	vcf1=${pid}-TuA_${pid}-NoA.var.vcf
	vcf2=${pid}-ReA_${pid}-NoA.var.vcf

	cnt1=$(grep -v ^# ${wd1}/${vcf1} | wc -l )
	cnt2=$(grep -v ^# ${wd1}/${vcf2} | wc -l )
	# vcf1=${pid}-TuA_${pid}-NoA_call_stats.txt.vcf.summary.exome_summary.csv.vcf.pass
	# vcf2=${pid}-ReA_${pid}-NoA_call_stats.txt.vcf.summary.exome_summary.csv.vcf.pass
	vcf1=${pid}-TuA_${pid}-NoA_call_stats.txt.vcf
	vcf2=${pid}-ReA_${pid}-NoA_call_stats.txt.vcf
	cnt3=$(grep -v ^# ${wd2}/${vcf1} | wc -l )
	cnt4=$(grep -v ^# ${wd2}/${vcf2} | wc -l )
	echo ${pid} ${cnt1} ${cnt2} ${cnt3} ${cnt4} >> temp_mutation_count
done < ${pidfile}

# get the number of mutated genes


	
##----------------------------
#step 2 randomize
# # 1. Total #bases in mRNA : 14,881,824,369
# 2. Total #bases in exons : 99,752,470
# 3. Total # bases in RefSeq Genes : 2,011,862,672




##----------------------------
#step 3 calculate probability 
Rscript ./get_Pvalue.r

