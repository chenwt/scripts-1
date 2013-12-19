#!/bin/bash
#By: J.He
#TODO: 

#mkfifo pipe1 pipe2 pipe3 
#grep "\.hg19\.seg" FILE_SAMPLE_MAP.txt > pipe1 &
#awk 'BEGIN{OFS="\t"}NR>1{print $1,$2".txt"}' pipe1 > pipe2 &
#cat pipe2 | cut -f2|sort | uniq -c |awk '$1>1{print $0}' 
#grep "\.hg19\.seg" FILE_SAMPLE_MAP_cnv_snparray.txt |grep -v "nocnv" >FILE_SAMPLE_MAP_cnv_snparray_hg19_cnvseg.txt
#
#awk -F"\t" 'NR>1{print $2}' FILE_SAMPLE_MAP_cnv_snparray_hg19_cnvseg.txt > pipe1 &
#awk -F"\t" 'NR>1{print $2}' FILE_SAMPLE_MAP_cnv_snparray_hg19_cnvseg.txt > pipe2 &
#awk -F"\t" 'NR>1{print $2}' FILE_SAMPLE_MAP_cnv_snparray_hg19_cnvseg.txt > pipe3 &
#cat pipe1 pipe2 pipe3 |sort |uniq -c |awk '$1>1{print $1}' |head 

#count=0
#for line in `ls *tar.gz`
#do
#   tar -xvzf $line
#   let count++
#   echo $count
#   mv FILE_SAMPLE_MAP.txt FILE_SAMPLE_MAP.txt$count
#done
#cat FILE_SAMPLE_MAP.txt* |sort|uniq > FILE_SAMPLE_MAP.txt

#grep "raw.copynumber.data.txt" FILE_SAMPLE_MAP.txt |wc -l 
#grep "birdseed.data.txt" FILE_SAMPLE_MAP.txt |wc -l 
#ls *raw.copynumber.data.txt |wc -l
#ls *birdseed.data.txt |wc -l

#grep ".birdseed.data.txt" FILE_SAMPLE_MAP_ALL.txt | awk 'BEGIN{
#      OFS="\t"
#    }
#    {if (length($2) > 25){
#      nf=split($2,a,","); print $1,a[nf]".txt"; }
#    else if (length($2) <=25 && length($2) >19)
#      print $1,$2".txt";
#    }' > /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/input_renameFiles.txt
##mv *birdseed.data.txt /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/
#wc -l /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/input_renameFiles.txt
#ls /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/*.birdseed.data.txt |wc -l
#
#grep "raw.copynumber.data.txt" FILE_SAMPLE_MAP_ALL.txt | awk 'BEGIN{
#      OFS="\t"
#    }
#    {if (length($2) > 25){
#      nf=split($2,a,","); print $1,a[nf]".txt"; }
#    else if (length($2) <=25 && length($2) >19)
#      print $1,$2".txt";
#    }' > input_renameFiles.txt
##
#wc -l input_renameFiles.txt
#ls *raw.copynumber.data.txt |wc -l 
#rm *ismpolish.data.txt
#rm *ismpolish.data.txt
#rm *byallele.copynumber.data.txt

#awk '{print substr($2,0,19)}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/input_renameFiles.txt > temp
#grep -v -f temp ../TCGA_barcode_all_in_cnv_meth_snp_EXP.txt > TCGA_barcode_to_download_again 
#rm temp
#mv FILE_SAMPLE_MAP.txt FILE_SAMPLE_MAP_ALL.txt
#tar -xzvf tcga_brca_cnv_level2_redownloaded.tar.tar.gz

#mv FILE_SAMPLE_MAP_ALL.txt FILE_SAMPLE_MAP_ALL.txt.temp
#cat FILE_SAMPLE_MAP.txt FILE_SAMPLE_MAP_ALL.txt.temp |sort|uniq > FILE_SAMPLE_MAP_ALL.txt
#wc FILE_SAMPLE_MAP_ALL.txt


~/scripts/projFocus/ceRNA/renameFiles.sh input_renameFiles.txt
#rm pipe1 pipe2 pipe3
