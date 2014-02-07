#!/bin/bash
#By: J.He
#TODO: 

###
# awk 'BEGIN{FS=OFS="\t"}NR>1&&$5!="NT"{print substr($2,0,19)}' brca_wgs_bam_summary_02042014.tsv > brca_wgs_bam_summary_02042014.tsv_TumorSample
# awk 'BEGIN{FS=OFS="\t"}NR>1&&$5=="NT"{print substr($2,0,19)}' brca_wgs_bam_summary_02042014.tsv > brca_wgs_bam_summary_02042014.tsv_NormalSample

# awk 'BEGIN{FS=OFS="\t"}NR>1&&$5!="NT"{print $2)}' brca_wgs_bam_summary_02042014.tsv > input_gtdownload_brca_wgs_bam_summary_02042014.tsv_TumorSample
grep -f brca_wgs_bam_summary_02042014.tsv_TumorSampleOverlapRNAseq  brca_wgs_bam_summary_02042014.tsv | awk -F"\t" '{print $2"\t"$17}' > input_gtBatch_v2.txt
min=0
for cnt in `seq 1 11` 
do
   max=
   min=`echo "$cnt*10"|bc` 
   max=
   awk -v n=$cntCut 'BEGIN{FS=OFS="\t"} NR<&&$5!="NT"{print $2)}' brca_wgs_bam_summary_02042014.tsv > input_gtdownload_brca_wgs_bam_summary_02042014.tsv


