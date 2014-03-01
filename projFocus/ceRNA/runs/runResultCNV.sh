#!/bin/bash
#! -cwd
#By: J.He
#Desp.: miscellic code for operations in this folder
#TODO: 

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/runs/runPipeline_step1.sh
# awk 'FNR!=1 ||NR==1' brca_cnvTumor_level2_naturePubNovalGene_02172014.mat brca_cnvTumor_level2_ucCeRNETCancerGene_02062014.mat |~/bin/sortxh |uniq > brca_cnvTumor_level2_combinedCG_02242014.mat 
# genelist='/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/tcgaPaper/tcga.16papers.genename' #the candidate cancer genes
# genelist='/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cancer.gene_UCceRNET.list' #the candidate cancer genes
# output=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_ucCeRNETCancerGene_02062014.mat
# getCNVMat $genelist $output
genelist='/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt'
output='/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_combinedCGDEG_Regulator_02282014.mat'
# getRegulatorCNVMat $genelist $output
# countSample $output
# python ~/scripts/projFocus/ceRNA/processData/getUniqCnvMat.py -i brca_cnvTumor_level2_combinedCGDEG_Regulator_02282014.mat -o brca_cnvTumor_level2_combinedCGDEG_Regulator_02282014_uniq.mat

