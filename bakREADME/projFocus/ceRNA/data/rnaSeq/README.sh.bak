#!/bin/bash
#By: J.He
#TODO: 

#for line in UNCID_1119153.TCGA-AN-A0AS-01A-11R-A00Z-07.110222_SN627_0063_B81FPVABXX.3_2.trimmed.annotated.gene.quantification.txt
#NCID_847434.TCGA-E2-A158-11A-22R-A12D-07.110629_UNC15-SN850_0090_AC03J3ABXX.8.trimmed.annotated.gene.quantification.txt ;do rdf UNCID_847434.TCGA-E2-A158-11A-22R-A12D-07.110629_UNC15-SN850_0090_AC03J3ABXX.8.trimmed.annotated.gene.quantification.txt >> input_getMat.txt; done

#for line in `cat input_getMat.txt `;do echo $line|awk -F"/" 'BEGIN{OFS="\t"}{split($NF,a,".");split(a[2],b,"-");print $0,$NF,b[3],substr(b[4],0,2)}' >>inpu_EdgeR_ColDesign.txt; done
#-------------------------------

#grep "gene.quantification.txt" FILE_SAMPLE_MAP.txt|awk 'BEGIN{OFS="\t"}{print $1,$2".txt"}' >input_renameFiles.txt
#head -2 input_renameFiles.txt
#awk '{print $2}' input_renameFiles.txt|sort |uniq -c |awk '$1>1{print $2}'|wc 
#~/scripts/projFocus/ceRNA/renameFiles.sh input_renameFiles.txt

#cut -f2 input_renameFiles.txt > input_makeMatRnaseql3.txt
#wc -l input_makeMatRnaseql3.txt 
#ls TCGA*txt |wc -l
#echo `date` >> log.README
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/makeMatRnaseql3.py -i input_makeMatRnaseql3.txt -o brca_rnaseq_rawCount_l3.mat >> log.README

#mv input_makeMatRnaseql3.txt input_makeMatRnaseql3.temp
#grep -f ../TCGA_barcode_all_in_cnv_meth_snp_EXP.txt input_makeMatRnaseql3.temp |sort > input_makeMatRnaseql3.txt 
echo `date` >> log.README
~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/makeMatRnaseql3.py -i input_makeMatRnaseql3.txt -o brca_rnaseq_rawCount_l3_geneSymbol.mat >>Rlog.README


