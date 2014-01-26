#!/bin/bash
#By: J.He
#$ -cwd
#TODO: 

### r script
##require(gdata)
##data = read.xls("GWAS_catalog_brca.xls")
##write.table(data,"GWAS_catalog_brca.txt",quote=F,row.names=F,sep="\t")
#colcount GWAS_catalog_brca.txt
#awk 'BEGIN{FS="\t"}{print $22}' GWAS_catalog_brca.txt|grep -e "^rs" |tr " ," "\n"|sed -e s/\s+//g |sort|uniq|sed 1d  >GWAS_catalog_brca_SNPid.txt 
#awk 'BEGIN{FS="\t"}{print $15}' GWAS_catalog_brca.txt|tr " \- " "\n" |sed 1d|sort|uniq|awk 'NR>11{print $0}' |tr ";" "\n"|sort|uniq > GWAS_catalog_brca_GeneName.txt

while read line
do
  grep -w $line ~/SCRATCH/database/projFocusRef/GenomeWideSNP_6.na22.annot.csv >> GWAS_catalog_brca_SNPid_Annot.csv
done < GWAS_catalog_brca_SNPid.txt 


