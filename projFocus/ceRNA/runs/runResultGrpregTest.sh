#!/bin/bash
#By: J.He
#TODO: 
#Desp: this is the running file for all coding testing in this folder

function qsub_Grplasso()
{
  chr=$1
  expfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno
 # snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr${chr}_KWtest.mat.anno.adjPass_1e-06.mat
  snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/test/brca_GWASCataLogGene_snp_KWtest.mat.anno.adjPass_1.0.mat
  cnvfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat
  cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh ${expfile} ${snpfile} ${cnvfile}" 
  $cmd
}







#~/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr21_KWtest.mat.anno.adjPass_1e-06.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat 
#~/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr21_KWtest.mat.anno.adjPass_1e-06.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat 

###get the table for genenumbers
#awk 'NR>1{print $2}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno | sort|uniq -c

##run with Dec 8th version
#qsub_Grplasso 21


for chr in 15 16 17 18 19
do
  if [ ! -d /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr ];
  then
    mkdir /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr 
  fi
  cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr
  echo `pwd` 
  qsub_Grplasso $chr
  echo $(date) 
  sleep 30m
done

#awk 'BEGIN{FS=":"}$NF>0.9&&FNR==1{print FILENAME,1-$NF}' grplasso_coeff_*.txt
#echo "#---------------"
#awk 'BEGIN{FS=":"}$NF<0.1&&FNR==1{print FILENAME,$NF}' grplasso_coeff_*.txt
