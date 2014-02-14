#!/bin/bash
#By: J.He
#TODO: 

cwd=`pwd`
cnvDir=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/CNV/report
cd $cnvDir
  grep PARGVC RN.xcnv.cnv.bed.geneAnnot |grep "^2" |cut -f5|sort|uniq  > $cwd/PARGVC_Re_chr2Dup_gene.txt

cd $cwd

grep -f MRs_02122014.txt PARGVC_Re_chr2Dup_gene.txt > PARGVC_Re_chr2Dup_gene_MRs_02122014.txt 
