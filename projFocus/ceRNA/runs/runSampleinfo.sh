#!/bin/bash
#By: J.He
#Desp:
#TODO: 
#$ -cwd

head -1 clinical_patient_brca.txt > TCGA_barcode_all_in_cnv_meth_snp_EXP_clinical.txt
while read line
do
  pid=`echo $line|awk '{print substr($1,0,12)}' `
  echo $pid
  grep -w $pid clinical_patient_brca.txt >> TCGA_barcode_all_in_cnv_meth_snp_EXP_clinical.txt
done < TCGA_barcode_all_in_cnv_meth_snp_EXP.txt


