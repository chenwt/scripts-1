#!/bin/bash
#By: J.He
#TODO: 

#mkfifo pipein pipeout
#cp input_renameFiles.txt input_renameFiles_old.txt
#grep "\.hg19\.seg\.txt" FILE_SAMPLE_MAP.txt | awk 'BEGIN{
#      OFS="\t"
#    }
#    {if (length($2) > 25){
#      nf=split($2,a,","); print $1,a[nf]".txt"; }
#    else if (length($2) <=25 && length($2) >19)
#      print $1,$2".txt";
#    }' > input_renameFiles.txt
#
#~/scripts/projFocus/ceRNA/renameFiles.sh input_renameFiles.txt
#awk '{split($2,a,"-");print a[4]}' input_renameFiles.txt |sort |uniq -c  
#head -1 input_renameFiles.txt 
###datat_start-------
#sample  Chromosome      Start   End     Num_Probes      Segment_Mean
#BEAUX_p_TCGA_b109_SNP_2N_GenomeWideSNP_6_A01_772082     1       61735   625470  34      0.0064
####data_end-------
#awk 'BEGIN{FS=OFS="\t"}{split($2,a,"-");print substr(a[4],0,2)}' input_renameFiles.txt | sort|uniq -c 
### get only tumor sample 
#---------

#awk 'BEGIN{FS=OFS="\t"}{
#    split($2,a,"-");
#    t=substr(a[4],0,2)
#    if(t~/^01/)
#      print a[3]
#    else
#      next
#    }' input_renameFiles.txt|sort|uniq >pid.tu.temp 
#
#wc -l pid.tu.temp
#awk 'BEGIN{FS=OFS="\t"}{
#    split($2,a,"-");
#    if(a[4]~/^1[01]/)
#      print a[3] 
#    else
#      next
#    }' input_renameFiles.txt|sort|uniq > pid.no.temp 
#wc -l pid.tu.temp
#grep -f pid.tu.temp pid.no.temp > pid.comm
#rm *temp
#wc -l pid.comm
#head pid.comm

#awk '{print $1"-01"}' pid.comm > temp.tu
#awk '{print $1"-1"}' pid.comm > temp.no
#ls TCGA*txt|tr "\t" "\n" |sort > temp
#rm input_getMat_tu.txt
#rm input_getMat_no.txt
#for i in  `cat pid.comm`
#do 
#  grep $i\-01 temp |head -1 >>input_getMat_tu.txt
#  egrep $i\-1[01] temp |head -1 >>input_getMat_no.txt
#done
#wc -l input_getMat_tu.txt
#wc -l input_getMat_no.txt
#

#head -450 input_getMat_tu.txt > input_getMat_tu_1.txt
#awk 'NR>450{print $0}' input_getMat_tu.txt > input_getMat_tu_2.txt
#for i in `cat input_getMat_tu_1.txt`
#do
#  wc -l $i
#done
~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/make_Mat_cnv_level3.py -i input_getMat_tu.txt -o brca_cnv_tumor.mat 



#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/compareFileRows.py -i brca_meth27k_matrix_l3.mat -m brca_meth450k_matrix_l3.mat -o brca_meth_level3_840.mat
#cat brca_meth_level3_840.mat.log
#wc -l brca_meth27k_matrix_l3.mat 
#awk 'NR==1{print NF}' brca_meth27k_matrix_l3.mat 
#wc -l brca_meth450k_matrix_l3.mat
#awk 'NR==1{print NF}' brca_meth450k_matrix_l3.mat
#wc -l brca_meth_level3_840.mat
#awk 'NR==1{print NF}' brca_meth_level3_840.mat
#rm pipein pipeout
