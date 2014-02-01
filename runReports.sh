#!/bin/bash
#By: J.He
#TODO: 

###----------------functions
##get variants count
awkWC(){
  awk 'END{print NR}' $1
}
getCnts(){
  dataSum=$1
  cntSNP=`awk 'END{print NR-1}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/brca_snp_tumor_731.mat.anno`
  cntCNV=`awk 'END{print NR-1}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/cnv/brca_gene_DEG_cnv_731.mat`
  cntSom=`awk 'END{print NR-1}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/somaticMutation/brca_som_selectedSample_level2.mat.anno `
  cntExp=`awk 'END{print NR-1}'  /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno `
  cntIndel=``
  echo -e "Type\tCount" > $dataSum
  echo -e "DEG\t$cntExp" >>$dataSum
  echo -e "SNP\t$cntSNP" >>$dataSum 
  echo -e "CNV\t$cntCNV" >>$dataSum
  echo -e "Somatic mutaiton\t$cntSom" >>$dataSum
  echo -e "Indel\t$cntIndel" >>$dataSum
  echo "#---DONE"
}
extractGene(){
  g=$1
  net=$2
  awk -v g=$g '$1==g||$2=g{print $0}' $net
}


##----------function---end


###--------RUN------
dataSum="data_summary.txt"
# getCnts $dataSum

##-----Gwas Catalog Gene run summary-----
##--Genes-
cntGene=`awkWC /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_allGeneName.txt`
cntf1=`/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_KWtest.mat.anno.adjPass_0.01.mat_GWASgene.mat`
cntf2=`/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat`
file="gcRUNSummary.txt"
echo -e "type\tcount" > $file
echo -e "GWAS catalog snp related Gene(GCgene)\t$cntGene" >> $file
echo -e "GCgene snp after KWtest(1e-2)\t$cntf1" >> $file
echo -e "GCgene snp after KWtest(1e-6)\t$cntf2" >> $file


