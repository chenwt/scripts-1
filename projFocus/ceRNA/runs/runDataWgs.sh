#!/bin/bash
#By: J.He
#TODO: 

###
# awk 'BEGIN{FS=OFS="\t"}NR>1&&$5!="NT"{print substr($2,0,19)}' brca_wgs_bam_summary_02042014.tsv > brca_wgs_bam_summary_02042014.tsv_TumorSample
# awk 'BEGIN{FS=OFS="\t"}NR>1&&$5=="NT"{print substr($2,0,19)}' brca_wgs_bam_summary_02042014.tsv > brca_wgs_bam_summary_02042014.tsv_NormalSample

# awk 'BEGIN{FS=OFS="\t"}NR>1&&$5!="NT"{print $2)}' brca_wgs_bam_summary_02042014.tsv > input_gtdownload_brca_wgs_bam_summary_02042014.tsv_TumorSample
# grep -f brca_wgs_bam_summary_02042014.tsv_TumorSampleOverlapRNAseq  brca_wgs_bam_summary_02042014.tsv | awk -F"\t" '{print $2"\t"$17}' > input_gtBatch_v2.txt
# min=0
# for cnt in `seq 1 11` 
# do
#    max=
#    min=`echo "$cnt*10"|bc` 
#    max=
#    awk -v n=$cntCut 'BEGIN{FS=OFS="\t"} NR<&&$5!="NT"{print $2)}' brca_wgs_bam_summary_02042014.tsv > input_gtdownload_brca_wgs_bam_summary_02042014.tsv
# done

# ls /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/TCGA-??-????-11*genome*csv.vcf > brca_wgs_normalbam_downloaded.txt
 # ls -1 ~/TCGA/BRCA/WGS/TCGA-??-????-11*bam >> brca_wgs_normalbam_downloaded.txt
# grep -f brca_wgs_normalbam_part2.txt brca_wgs_bam_summary_02042014.tsv| awk -F"\t" '{print $2"\t"$17}' > input_gtBatch_normalBam_part2.txt

# grep -w  HMS-RK ~/SCRATCH/projFocus/ceRNA/data/wgs/brca_wgs_bam_summary_02042014.tsv |cut -f2 > temp
# grep -f temp input_gtBatch_v2.txt >> input_gtBatch_v2_broadReference.txt

grep -v -w  HMS-RK ~/SCRATCH/projFocus/ceRNA/data/wgs/brca_wgs_bam_summary_02042014.tsv |cut -f2 > temp
grep -f temp input_gtBatch_v2.txt >> input_gtBatch_v2_WU.txt

rm temp
