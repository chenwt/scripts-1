#!/bin/bash
#$ -cwd
#By: J.He
#TODO: 
#Desp: this is the running file for all coding testing in this folder
##run on selected know BRCA genes

###---prepare for data

function extractOneData4Gene()
{
  #function to get record for each gene giving one data
  #called by extractAllData4Gene
  gene=$1
  file=$2
  outfile=`echo $file|awk 'BEGIN{FS="/"}{print $NF}'`
  outfile=$gene"_"$outfile
  awk -v gene=$gene 'NR==1||$1==gene{print $0}' ${file} > $outfile
}

function extractAllData4Gene()
{
  #given chr and gene names, get expression, snp, cnv, and somatic information, output to current working DIR
  chr=$1
  gene=$2
  expfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno
  snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr${chr}_KWtest.mat.anno.adjPass_1e-06.mat
  cnvfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat
  somfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/som/brca_somForDeg.mat
  extractOneData4Gene $gene ${expfile} 
  extractOneData4Gene $gene ${snpfile} 
  extractOneData4Gene $gene ${cnvfile} 
  extractOneData4Gene $gene ${somfile}
}

function pSystime()
{
  #function to print current system time
 t=$(date) 
  echo -e "System time:\t"$t
}
##-----------------------------------------------------------
##-------initiation
##--global variables
rootwd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/fTest/test
snpDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp
somDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/som
cnvDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv
expDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp
genelist=$rootwd/knowledgeBase/knowGene.txt
outputFile="fTest_pval.txt" ##filenames cannot be changed since it's in the R script

##--preprocessing
##replace knowGene.txt with any gene list in future run
if [ ! -f $wd/knowBrcaGene_DEG.anno.txt ] 
then
 grep -wf $genelist $expDir/brca_exp_l3_731_DEG.mat.signleTSS.anno.GeneName > $wd/knowBrcaGene_DEG.txt
 grep -wf $wd/knowBrcaGene_DEG.txt ~/SCRATCH/database/projFocusRef/gene_annotation_hg19_unique_start_end.txt >$wd/knowBrcaGene_DEG.anno.txt
 rm  $wd/knowBrcaGene_DEG.txt
fi 


###------run ftest.r for input genelist
function runFtest() {
    while read line 
    do
      chr=`echo $line|awk '{print $2}'`
      gene=`echo $line|awk '{print $1}'`
      echo -e "chr:\t"$chr"\tgene:\t"$gene
      echo -e "get files...."
      pSystime 
      extractAllData4Gene $chr $gene  
      echo -e "run r...."
      pSystime
      cntSnp=`awk 'END{print NR}' ${gene}_brca_gene_snp_chr${chr}_KWtest.mat.anno.adjPass_1e-06.mat`
      cntCnv=`awk 'END{print NR}' ${gene}_brca_somForDeg.mat `
      cntSom=`awk 'END{print NR}' ${gene}_brca_gene_DEG_cnv_731.mat `
      flag=1
      if [ $cntSnp -eq 1 ] ; then  exit; fi
      if [ $cntCnv -gt 1 ] & [ $cntSom -gt 1 ] ; then
          /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript ~/scripts/projFocus/ceRNA/ftest.r --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_brca_gene_snp_chr${chr}_KWtest.mat.anno.adjPass_1e-06.mat --som ${gene}_brca_somForDeg.mat --cnv ${gene}_brca_gene_DEG_cnv_731.mat --gene $gene 
      else
        echo -e "$gene\tNA\tNA" >> $outputFile
      fi
      pSystime 
      echo -e "cleaning file..."
      pSystime 
      rm ${gene}*
      
    done < $wd/knowBrcaGene_DEG.anno.txt
}
#runFtest

##---get significant genes(snp significantly contribute)
cd $wd
awk '$2<0.1{split($1,a,"_");print a[1]}' $outputFile > $outputFile.sig

###----run fTest.r for all input genes
#~/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr21_KWtest.mat.anno.adjPass_1e-06.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat 
#~/scripts/projFocus/ceRNA/qsub_grpLassoSNP_v1.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr21_KWtest.mat.anno.adjPass_1e-06.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat 

###get the table for genenumbers
#awk 'NR>1{print $2}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno | sort|uniq -c


##run with version3 grplasso.R
#for chr in 15 16 17 18 
#do
#  if [ ! -d /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr ];
#  then
#    mkdir /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr 
#  fi
#  cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr$chr
#  echo `pwd` 
#  qsub_ridgeFtest $chr
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
#  qsub_ridgeFtest $chr
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
#  qsub_ridgeFtest $chr
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

#grpCNV 21
#awk 'BEGIN{FS=":"}$NF>0.9&&FNR==1{print FILENAME,1-$NF}' chr*/grplasso_coeff_*.txt
#echo "#---------------"
#awk 'BEGIN{FS=":"}$NF<0.1&&FNR==1{print FILENAME,$NF}' chr*/grplasso_coeff_*.txt |less

##number of job submitted
#wc chr*/log/qsub*

