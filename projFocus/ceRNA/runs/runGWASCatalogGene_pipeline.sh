#!/bin/bash
#By: J.He
#TODO: 
#$ -cwd

###---rerun for know GWAS genes and GWAS snps
####------generateing matrix files

##step 1
#cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/makeMat_SNP.py -i input_makMat_SNP_tumorComm.txt -o brca_snp_tumor >> log.RUN
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/annotSNP_v2.py -i brca_snp_tumor_731.mat -d ~/SCRATCH/database/projFocusRef/annot_GenomeWideSNP_6_hg19_5col.txt -o brca_snp_tumor_731.mat.anno &

####---------END-----------------------------------


###-----pre-processing files-----------------------
##step 1
#annotate differential expression gene expression--ran long time ago
#/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/annotSNP_v2.py -i brca_rnaseq_rawCount_l3_geneSymbol_DEG_Mat.txt -d /ifs/scratch/c2b2/ac_lab/jh3283/database/projFocusRef/refseq_gene_tss_geneSymbol.signleTSS.5cols -o brca_exp_l3_731_DEG.mat.singleTSS.anno


##step 2-
#cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/test
##extract knowSNPs
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/extractSNPmat.py -i ~/SCRATCH/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_SNPid.txt -m /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/brca_snp_tumor_731.mat.anno -o brca_snp_tumor_731_GWASCatalogSNP.mat.anno 

##step 3.1---prepare snp file KW test and annotated with GeneName: have to run this step to associate SNP with gene by TSS position
#cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/test
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731_GWASCatalogSNP.mat.anno  -e  ../brca_exp_l3_731_DEG.mat.singleTSS.anno -o brca_GWASCataLogGene_snp_KWtest.mat.anno -j 1
#
#cat brca_GWASCataLogGene_snp_KWtest.mat.anno*log
#rm brca_GWASCataLogGene_snp_KWtest.mat.anno*log
#rm brca_GWASCataLogGene_snp_KWtest.mat.anno*snp


##step 3.2 
##prepre expression file
 

##step 3.3-
##prepare cnv file
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/mergeExpCNV.py -c input_for_mergeSMC_v1_cnv.txt.sorted  -e /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno  -o brca_gene_DEG_cnv_731.mat

##step 3.4
##prepare som file
#cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/som
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/mergeExpSomMutation.py -s /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/somaticMutation/brc


##step 6
##run Ftest on each gene
# cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/fTest/test
# awk 'NR>1{print $1}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno.GWASCataGene.mat.anno |sort |uniq> gene.list

#--step 7
##run grp lasso model
# cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grplasso/test




