#!/bin/bash
#By: J.He
#TODO: 
##Desp: after gatk annoatation, delet file in TCGA folder

resultDir=$1
dataDir=$2
cnt=0
for vcf in `ls $resultDir/*_dindel_ouput.variantCalls.VCF`
do
  bamName=`echo $vcf |awk -F/ '{gsub("_dindel_ouput.variantCalls.VCF","",$NF);print $NF}'`
  # echo $bamName
  if [ -f $dataDir/$bamName ]; then
    echo -e "Delete\t $dataDir/$bamName"
    rm $dataDir/$bamName
    rm $dataDir/$bamName.bai
    ((cnt+=1))
  fi
done

echo -e "$cnt BAM file deleted"
    

