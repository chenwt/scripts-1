#!/bin/bash
#By: J.He
#TODO: 

mkfifo pipein pipeout
#awk '$2!~/^\s+/{print substr($2,0,19)}' FILE_SAMPLE_MAP_meth.txt |sort |head 
#head -2 broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt |tail -1 > pipein &
#awk -F"\t" '{if (NR==1) { for (i=1; i<=NF; i++) {printf i; for (j=1; j<=int(length($i)/8); j++) {printf "\t"}; if (i==NF) {printf "\n"} else {printf "\t"}};print $0 } else {print $0}}' pipein | less -S

#awk -F"\t" 'BEGIN{OFS="\t"}$31~/[A-Za-z1-9]/{print $31,$2}' broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt  > FILE_MAP_BRCA.Genome_Wide_SNP_6_birdseed.txt
#awk -F"\t" 'BEGIN{OFS="\t"}$76~/[A-Za-z1-9]/{print $76,$2}' broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt  > FILE_MAP_BRCA.Genome_Wide_SNP_6_cnv_seg_hg19.txt
#awk -F"\t" 'BEGIN{OFS="\t"}$40~/[A-Za-z1-9]/{print $40,$2}' broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt  > FILE_MAP_BRCA.Genome_Wide_SNP_6_copynumber_raw.txt
#
##grep "\.hg19\.seg" FILE_SAMPLE_MAP_cnv_snparray.txt |grep -v "nocnv" >FILE_SAMPLE_MAP_cnv_snparray_hg19_cnvseg.txt
#
#mkfifo pipe1 pipe2 pipe3
#awk -F"\t" 'NR>1{print substr($2,0,19)}' FILE_MAP_BRCA.Genome_Wide_SNP_6_cnv_seg_hg19.txt|sort|uniq > pipe1 &
#awk -F"\t" 'NR>1&&$1!~/^\s+/{print substr($2,0,19)}' FILE_SAMPLE_MAP_meth.txt |sort|uniq > pipe2 &
#awk -F"\t" 'NR>1{print substr($2,0,19)}' FILE_MAP_RCA.Genome_Wide_SNP_6_birdseed.txt |sort|uniq > pipe3 &
###head -2 pipe1
###echo -e "-------"
###head -2 pipe2
###echo -e "-------"
###head -2 pipe3
###echo -e "-------"
#cat pipe1 pipe2 pipe3 |sort |uniq -c |awk '$1==3{print $2}' > TCGA_barcode_all_in_cnv_meth_snp.txt
#rm pipe1 pipe2 pipe3

##cp rnaSeq/FILE_SAMPLE_MAP.txt FILE_SAMPLE_MAP_rnaseq.txt
grep "gene\.quantification\.txt" FILE_SAMPLE_MAP_rnaseq.txt |awk -F"\t" '{print substr($2,0,19)}' > pipein &
cat TCGA_barcode_all_in_cnv_meth_snp.txt pipein |sort|uniq -c |awk '$1==2{print $2}' > TCGA_barcode_all_in_cnv_meth_snp_EXP.txt
wc TCGA_barcode_all_in_cnv_meth_snp_EXP.txt
#awk '{print substr($1,14,3)}' TCGA_barcode_all_in_cnv_meth_snp_EXP.txt|sort|uniq -c 
#cat pipein |awk '{print substr($1,14,3)}'|sort|uniq -c 
#cat pipein |awk '{print substr($1,14,3)}'|sort|uniq -c 
#rm TCGA_barcode_all_in_cnv_meth_snp_EXP.txt



rm pipein pipeout

