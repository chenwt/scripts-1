#!/bin/bash
#By: J.He
#TODO: 


#ln -s ../linkSNP_EXP/edgeR_brca_l3DEGMat.txt.anno.comm brca_exp_edgeR_DEG.mat.anno.comm
#ln -s ../linkSNP_EXP/brca_snp_level3_839.mat.anno.comm brca_snp_level3_839.mat.anno.comm
#ln -s ../meth/brca_meth_matrix_l3.mat.anno.Mval.comm brca_meth_matrix_l3.mat.anno.Mval.comm

### generate test file for grplasso

#awk '$2==22{print $0}' brca_exp_edgeR_DEG.mat.anno.comm        >/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/test_exp_chr22.mat
#awk '$2==22{print $0}' brca_snp_level3_839.mat.anno.comm       >/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/test_snp_chr22.mat
#awk '$2==22{print $0}' brca_meth_matrix_l3.mat.anno.Mval.comm  >/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/test_meth_chr22.mat

####-------prepare_for_the_input_file------
#make sure the order of all columns are same! 
#head -1 brca_exp_edgeR_DEG.mat.anno.comm|tr "\t" "\n" > header_exp
#head -1 brca_snp_level3_839.mat.anno.comm |tr "\t" "\n" > header_snp
#head -1 brca_meth_matrix_l3.mat.anno.Mval.comm |tr "\t" "\n" > header_meth
#cat cnv_level3/input_getMat_tu.txt cnv_level3/input_getMat_no.txt  > header_cnv
normalIndex (){
awk 'NR>4{old=$1;split($1,a,".");
   if(a[4]~/0[13]/)   
     print $1
   else
     next
    }' $1 > $1.tu 
wc $1.2cols
cut -f1 $1.2cols| cut -d. -f4 |sort |uniq -c
}
#normalIndex header_exp 
#normalIndex header_snp
#normalIndex header_meth

###--------------perpare_for_the_cnv.mat_file
#awk '{gsub("-",".",$1);print substr($1,0,19)}' header_cnv > temp.header_cnv
#awk '{print substr($1,0,19)}' header_exp.2cols > temp.header_exp
#grep -f temp.header_exp temp.header_cnv > final.sample.tu 
#grep -f final.sample.tu header_exp  > final.header_exp
#grep -f final.sample.tu header_snp  > final.header_snp
#grep -f final.sample.tu header_meth > final.header_meth
#sed "s/\./-/g" final.sample.tu > final.sample.tu.temp

#head -1 final.header_exp |tr "\t" "\n"|awk '{print substr($1,0,19)}'|sort > temp1
#head -1 final.header_snp |tr "\t" "\n" |awk '{print substr($1,0,19)}'|sort > temp2
#head -1 final.header_meth |tr "\t" "\n" |awk '{print substr($1,0,19)}'|sort > temp3
#comm -3 temp1 temp2
#comm -3 temp2 temp3
#comm -3 temp1 temp3


#head -5 final.header_*
#wc final.header_*
#cp final.header_cnv input_for_mergeSMC_v1_cnv.txt
#for line in `ls final.header_*`
#do
#  sort $line > $line.sorted
#done
#
#if [ -f final.header_cnv ]
#then
#  rm final.header_cnv
#fi
#awk '{print substr($1,0,19)}' final.header_exp.sorted > temp.final.header_exp.sorted
#for pid in `cat temp.final.header_exp.sorted`
#do
#  grep $pid header_cnv >> final.header_cnv
#done
#mv final.header_cnv input_for_mergeSMC_v1_cnv.txt.sorted
#
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/getCols.py -i brca_exp_edgeR_DEG.mat.anno.comm   -c final.header_exp.sorted -o brca_exp.mat.anno.comm.finalTuSample 
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/getCols.py -i brca_snp_level3_839.mat.anno.comm  -c final.header_snp.sorted -o brca_snp.mat.anno.comm.finalTuSample  
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/getCols.py -i brca_meth_matrix_l3.mat.anno.Mval.comm   -c final.header_meth.sorted -o brca_meth.mat.anno.comm.finalTuSample 
##double checking....
#head -1 brca_exp.mat.anno.comm.finalTuSample |tr "\t" "\n"|awk '{print substr($1,0,19)}'|sort > temp1
#head -1 brca_snp.mat.anno.comm.finalTuSample |tr "\t" "\n" |awk '{print substr($1,0,19)}'|sort > temp2
#head -1 brca_meth.mat.anno.comm.finalTuSample |tr "\t" "\n" |awk '{print substr($1,0,19)}'|sort > temp3
#comm -3 temp1 temp2
#comm -3 temp2 temp3
#comm -3 temp1 temp3

##-----------------prepare_for_gene_CNV_mat-------
### create softlink of cnv level3 data
for file in `cat input_for_mergeSMC_v1_cnv.txt.sorted`
do
  if [ ! -f $file ] 
  then
    ln -s cnv_level3/$file $file
  fi
done 
~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/mergeExpCNV.py -c input_for_mergeSMC_v1_cnv.txt.sorted  -e /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno  -o brca_gene_DEG_cnv_731.mat
find . -type l -print |grep TCGA|xargs rm

#####------------run_main_code-------
##~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/getSingleTSS_exp.py brca_exp.mat.anno.comm.finalTuSample brca_exp.mat.anno.comm.finalTuSample.singleTSS
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/merge_snp_meth_cnvl3_v1.py -s brca_snp.mat.anno.comm.finalTuSample -m brca_meth.mat.anno.comm.finalTuSample -c input_for_mergeSMC_v1_cnv.txt.sorted -g brca_exp.mat.anno.comm.finalTuSample.singleTSS -o brca_gene_snp_meth_cnv.mat.finalTuSample

#/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/merge_snp_meth_cnvl3_v1.py




###prepare_for__regression___test
#rm knownBRCA.gene.entrezid
#for gene in `cat knownBRCA.gene`
#do
#  echo $gene
#  grep "$gene"  gene_RefSeqGene |head -1 
#  grep "$gene"  gene_RefSeqGene |head -1 | awk '{print $2}' >> knownBRCA.gene.entrezid
#done

#awk -f knownBRCA.gene.entrezid brca_gene_snp_meth_cnv.mat.finalTuSample > test_knowBRCA.gene_snp_meth_cnv.mat.finalTuSample
#awk 'FNR==NR {a[$1];next} ($1 in a)' knownBRCA.gene.entrezid brca_gene_snp_meth_cnv.mat.finalTuSample > test_knowBRCA.gene_snp_meth_cnv
#awk 'FNR==NR {a[$1];next} ($1 in a)' knownBRCA.gene.entrezid brca_exp.mat.anno.comm.finalTuSample.singleTSS  > test_knowBRCA.gene_exp
#readlink -f test_knowBRCA.gene_exp

#Rscript  /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/grpreg.r /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test_knowBRCA.gene_snp_meth_cnv
