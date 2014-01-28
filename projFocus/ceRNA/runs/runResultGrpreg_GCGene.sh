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

function extractAllData4Gene_GCgene()
{
  #only changes the snp data searching space
  gene=$1
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

function runFtest() {
    while read line 
    do
      gene=`echo $line|awk '{print $1}'`
      echo -e "gene:\t"$gene
      #echo -e "get files...."
      #pSystime 
      extractAllData4Gene_GCgene $gene  
      #pSystime
      cntSnp=`awk 'END{print NR}' ${gene}_brca_GWASCataLogGene_snp_KWtest.mat.anno.adjPass_1.0.mat`
      if [ $cntSnp -eq 1 ] ; then 
	echo -e "$gene" >> ${outputFile}_nosnp.txt
       	continue 
      fi
      cntSom=`awk 'END{print NR}' ${gene}_brca_somForDeg.mat `
      cntCnv=`awk 'END{print NR}' ${gene}_brca_gene_DEG_cnv_731.mat `
      ##neet to add cntIndel in future

      echo -e "cntSnp \t"${cntSnp}"\tcntSom:\t"${cntSom}"\tcntCnv:\t"$cntCnv
      
      if [[ $cntCnv -gt 1 && $cntSom -eq 1 ]] ; then
      ##case 2 snp + cnv
      echo -e "snp cnv...."
          /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript ~/scripts/projFocus/ceRNA/ftest_v2.r --type 2 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_${snpFile} --som ${gene}_brca_somForDeg.mat --cnv ${gene}_brca_gene_DEG_cnv_731.mat --gene $gene --out ${outputFile}_snp_cnv.txt 

      elif [[ $cntCnt -eq 1 && $cntSom -gt 1 ]] ; then
      ##case 3 snp + som
      echo -e "snp som...."
          /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript ~/scripts/projFocus/ceRNA/ftest_v2.r --type 3 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_${snpFile} --som ${gene}_brca_somForDeg.mat  --out ${outputFile}_snp_som.txt 
          	
      elif [[ $cntCnv -gt 1 && $cntSom -gt 1 ]] ; then
      ##case 1 snp + cnv + som
      echo -e "snp cnv som...."
          /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript ~/scripts/projFocus/ceRNA/ftest_v2.r --type 1 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_${snpFile} --som ${gene}_brca_somForDeg.mat --cnv ${gene}_brca_gene_DEG_cnv_731.mat --gene $gene --out ${outputFile}_snp_cnv_som.txt 

      else  
	echo "No cnv, som... for $gene"
        echo -e "$gene" >> ${outputFile}_snp.txt
      fi
      rm ${gene}*
      
    done < $1

    for file in `ls ${outputFile}_snp_*txt`
    do
        awk '$2<0.1{split($1,a,"_");print a[1]"\t"$2}' $file|sort -k 2,2n|uniq > $file.sig
    done

    echo "#-----END------"
}

##--global variables
expfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno
snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/test/brca_GWASCataLogGene_snp_KWtest.mat.anno.adjPass_1.0.mat
cnvfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat
somfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/som/brca_somForDeg.mat
#cut -f1 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/fTest/test/fTest_pval_snp_cnv_som.txt.sig  > gene.list
#while read gene
#do
# extractAllData4Gene_GCgene  ${gene}
#done < gene.list
extractAllData4Gene_GCgene 'ESR1'
###------initiation
#rootwd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/
#wd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test
#genelist=$wd/gene.list
#outputFile="grplasso_coeff" 
#snpFile=brca_GWASCataLogGene_snp_KWtest.mat.anno.adjPass_1.0.mat
###------run ftest.r for input genelist

#runFtest $genelist 
#
###---get significant genes(snp significantly contribute)
#function grpCNV (){
#  cwd=`pwd`
#  for gene in `ls $cwd/chr$1/grplasso_*RData |awk 'BEGIN{FS="_"}{split($NF,a,"\.");print a[1]}'`
#  do
#    echo -e "#------chr\t$1\tgene\t"$gene
#    cd $cwd/chr$1
#    Rscript ~/scripts/projFocus/ceRNA/grpLassoSNP_v3_cnv.r temp/input_snp_$gene temp/input_exp_$gene temp/input_cnv_$gene 0
#  done
#}
#
