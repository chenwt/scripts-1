#!/bin/bash
#By: J.He
#TODO: 
#Desp: this is the running file for all coding testing in this folder

#~/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr21_KWtest.mat.anno.adjPass_1e-06.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat 
#~/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr21_KWtest.mat.anno.adjPass_1e-06.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat 

###get the table for genenumbers
#awk 'NR>1{print $2}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno | sort|uniq -c

function qsub_Grplasso()
{
  chr=$1
  expfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno
  snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr${chr}_KWtest.mat.anno.adjPass_1e-06.mat
  cnvfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat
  cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh ${expfile} ${snpfile} ${cnvfile}" 
  $cmd
}
##run with cnv data,keep the original directory structure
function qsub_GrplassoCnv()
{
  chr=$1
  expfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno
  snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr${chr}_KWtest.mat.anno.adjPass_1e-06.mat
  cnvfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat
  cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v3.sh ${expfile} ${snpfile} ${cnvfile}" 
  $cmd
}


##run with Dec 8th version
#qsub_Grplasso 21

##run with version3 grplasso.R
#for chr in 15 16 17 18 
#do
#  if [ ! -d /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr ];
#  then
#    mkdir /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr 
#  fi
#  cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr
#  echo `pwd` 
#  qsub_Grplasso $chr
#  echo $(date) 
#  sleep 30m
#done

#run v3 without plot for chromsom 1-3, 4-7, 8-10,11-14,19
#for chr in `seq 20 22` X Y  
#do
#  if [ ! -d /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr ];
#  then
#    mkdir /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr 
#  fi
#  cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr
#  #when rerun 20:22 X Y
#  #-----
#  if [! -d resultV2 ];
#  then 
#    mkdir resultV2
#    mv grplasso_coeff_* resultV2/
#  fi
#  #-----
#  echo `pwd` 
#  qsub_Grplasso $chr
#  echo $(date) 
#  sleep 20m
#done

##run v3_cnv version with qsub_v3 for all chromosomes 
#test on Y first
#for chr in Y  
#do
#  if [ ! -d /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr ];
#  then
#    mkdir /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr 
#  fi
#  cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr
#  echo `pwd` 
#  qsub_Grplasso $chr
#  echo $(date) 
#  #sleep 20m
#done

function grpCNV (){
  cwd=`pwd`
  for gene in `ls $cwd/chr$1/grplasso_*RData |awk 'BEGIN{FS="_"}{split($NF,a,"\.");print a[1]}'`
  do
    echo -e "#------chr\t$1\tgene\t"$gene
    cd $cwd/chr$1
    Rscript ~/scripts/projFocus/ceRNA/grpLassoSNP_v3_cnv.r temp/input_snp_$gene temp/input_exp_$gene temp/input_cnv_$gene 0
  done
}

grpCNV 21
#awk 'BEGIN{FS=":"}$NF>0.9&&FNR==1{print FILENAME,1-$NF}' chr*/grplasso_coeff_*.txt
#echo "#---------------"
#awk 'BEGIN{FS=":"}$NF<0.1&&FNR==1{print FILENAME,$NF}' chr*/grplasso_coeff_*.txt |less

##number of job submitted
#wc chr*/log/qsub*

