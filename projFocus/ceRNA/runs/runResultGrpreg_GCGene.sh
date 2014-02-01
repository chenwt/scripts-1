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
 ##--global variables
  expfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno
  #snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_KWtest.mat.anno.adjPass_0.01.mat_GWASgene.mat
  snpfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat
  cnvfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat
  somfile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/som/brca_somForDeg.mat

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
function grpLasso (){
###-------under development
  cwd=`pwd`
  #for gene in `ls $cwd/chr$1/grplasso_*RData |awk 'BEGIN{FS="_"}{split($NF,a,"\.");print a[1]}'`
  for gene in ``
  do
    echo -e "#------chr\t$1\tgene\t"$gene
    cd $cwd/chr$1
   /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r --exp ESR1_brca_exp_l3_731_DEG.mat.singleTSS.anno --snp input_test_reg_snp.mat --cnv ESR1_brca_gene_DEG_cnv_731.mat --som ESR1_brca_somForDeg.mat --type 1 --plot 0 --out grplasso_
  done
}


function runGrplasso() {
    while read line 
    do
      gene=`echo $line|awk '{print $1}'`
      echo -e "gene:\t"$gene
      #echo -e "get files...."
      #pSystime 
      extractAllData4Gene_GCgene $gene  
      #pSystime
      cntSnp=`awk 'END{print NR}' ${gene}_${snpFile}`
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
       cmd="/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r  --plot 0 --type 2 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_${snpFile}  --cnv ${gene}_brca_gene_DEG_cnv_731.mat --gene $gene --out ${outputFile}"

      elif [[ $cntCnv -eq 1 && $cntSom -gt 1 ]] ; then
      ##case 3 snp + som
      echo -e "snp som...."
         cmd="/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r  --plot 0 --type 3 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_${snpFile} --som ${gene}_brca_somForDeg.mat  --gene $gene --out ${outputFile}"
          	
      elif [[ $cntCnv -gt 1 && $cntSom -gt 1 ]] ; then
      ##case 1 snp + cnv + som
	echo -e "snp cnv som...."
         cmd="/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r  --plot 0 --type 1 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno  --snp ${gene}_${snpFile} --som ${gene}_brca_somForDeg.mat --cnv ${gene}_brca_gene_DEG_cnv_731.mat --gene $gene --out ${outputFile}"
      elif [[ $cntCnv -eq 1 && $cntSom -eq 1 ]] ; then
	##case 4 snp.
         echo -e "snp ...."
         cmd="/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r  --plot 0 --type 4 --exp ${gene}_brca_exp_l3_731_DEG.mat.singleTSS.anno --snp ${gene}_${snpFile} --out ${outputFile}"

      else  
	echo "No cnv, som, snp for $gene"
	echo " ERROR---"
        echo -e "$gene" >> ${outputFile}.noDataGene
      fi
      $cmd &
      #echo $cmd|qsub -l mem=8g,time=20:: -N grp_$gene -e ./log -o ./log -cwd > qsub.log
      # tail -1 qsub.log
      
    done < $1
    echo "#-----END------"
}

function getGenelist(){
  awk '{print $1}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/fTest/test/fTest_pval_all__snp_cnv_som.txt.sig > gene.list 
}

##filtering and save rda file
function runFilterGrpLasso(){
  cutoff=$1
  /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterGrplasso.r --file grplasso_coeff --cut $cutoff --out gcGenes_GeneVarNet_${cutoff}
}

###-----------------
####-----------function---end--------------------------------

##------initiation
rootwd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test
#snpFile=brca_gene_snp_KWtest.mat.anno.adjPass_0.01.mat_GWASgene.mat
# ln -s /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat 
snpFile=brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat

getGenelist 
while read gene
do
  extractAllData4Gene_GCgene ${gene}
done < gene.list

genelist=$wd/gene.list
outputFile="grplasso_coeff" 
runGrplasso $genelist 
runFilterGrpLasso 0.05

rm ${gene}*
##---get significant genes(snp significantly contribute)
 #/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r --exp ESR1_brca_exp_l3_731_DEG.mat.singleTSS.anno --snp input_test_reg_snp.mat --cnv ESR1_brca_gene_DEG_cnv_731.mat --som ESR1_brca_somForDeg.mat --type 1 --plot 0 --out grplasso_

# mkdir run001
# mv grplasso* run001 
