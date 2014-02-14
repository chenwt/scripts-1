#!/bin/bash
#By: J.He
#TODO: 
##Desp: after gatk annoatation, delet file in TCGA folder

resultDir=$1
dataDir=$2
cnt=0
for vcf in `ls $resultDir/*var.vcf.gatk.vcf`
do
  bamName=`echo $vcf |awk -F/ '{gsub(".var.vcf.gatk.vcf","",$NF);print $NF}'`
  echo $bamName
  if [ -f $dataDir/$bamName ]; then
    echo -e "Delete\t $dataDir/$bamName"
    rm $dataDir/$bamName
    ((cnt+=1))
  fi
done

echo -e "$cnt BAM file deleted"
    

