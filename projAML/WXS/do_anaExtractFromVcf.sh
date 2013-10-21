#!/bin/bash
#By: J.He
#TODO: modify path for muTect result 

sd=~scripts/projAML/
#wd=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/callVars/result/
#od=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/report/
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/muTect/varAccepted/
od=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/muTect/report/
pidf=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/PID_16.txt

cd ${wd}

while read  pid
do
    #awk '!/^#/ {print $1,$2}' ${pid}-ReA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv.vcf > RN.temp
    #awk '!/^#/ {print $1,$2}' ${pid}-TuA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv.vcf > TN.temp
    #comm -13 <(awk 'BEGIN{OFS=FS="\t"} !/^#/ {print $1,$2}' ${wd}${pid}-TuA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv.vcf | sort) <( awk 'BEGIN{OFS=FS="\t"} !/^#/ {print $1,$2}' ${wd}${pid}-ReA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv.vcf| sort) > RN_minus_TN.temp 
   # head -1 ${wd}${pid}-ReA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv > ${od}${pid}_RN-TN.txt
    #grep -f RN_minus_TN.temp ${wd}${pid}-ReA_${pid}-NoA.var.vcf.somaticfiltered.vcf.summary.exome_summary.csv > ${od}${pid}_RN-TN.txt

    awk '!/^#/ {print $1,$2}' ${pid}-ReA_${pid}-NoA_call_stats.txt.vcf.summary.genome_summary.csv.vcf.pass  > RN.temp
    awk '!/^#/ {print $1,$2}' ${pid}-TuA_${pid}-NoA_call_stats.txt.vcf.summary.genome_summary.csv.vcf.pass > TN.temp
    comm -13 <(awk 'BEGIN{OFS=FS="\t"} !/^#/ {print $1,$2}' ${wd}${pid}-TuA_${pid}-NoA_call_stats.txt.vcf.summary.genome_summary.csv.vcf.pass | sort) <( awk 'BEGIN{OFS=FS="\t"} !/^#/ {print $1,$2}' ${wd}${pid}-ReA_${pid}-NoA_call_stats.txt.vcf.summary.genome_summary.csv.vcf.pass | sort) > RN_minus_TN.temp 
    grep -f RN_minus_TN.temp ${wd}${pid}-ReA_${pid}-NoA_call_stats.txt.vcf.summary.genome_summary.csv.pass > ${od}${pid}_RN-TN.txt

done < ${pidf}
    #ls *NoA*\.vcf|grep -f ${pidf} > vcf2procee.temp




