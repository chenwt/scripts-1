#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 


#-----------function
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/runs/runPipeline_step1.sh 
source ~/scripts/projFocus/ceRNA/model/func4DoModel.sh  

srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/
RSCRIPT=/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript
PYTHON=~/tools/python/Python_current/python
#-----------execute
END="#---[THE END]---"
#------------step1-prepareData

##-----------# #----------##-----------# #-----------sec
##-----------# #----------##-----------# 
###---second run for nature 12912 genelist

genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list
cnvmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_naturePubNovalGene_02172014.mat
# if [ ! -f $cnvmat ]; then
#   getCNVMat $genelist $cnvmat 
# fi
#---
tumorMethDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/methlevel3/
normalMethDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/methlevel3Normal/
genelistAnno=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list.geneSingleStartEnd
outMethTumor=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth/brca_methTumor_nature12912NovelGene_02172014.mat
outMethNormal=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth/brca_methNormal_nature12912NovelGene_02172014.mat

# if [ ! -f $outMethTumor ]; then
#   genMethMat $tumorMethDir $genelistAnno $outMethTumor  
# fi
# if [ ! -f $outMethNormal ]; then
#   genMethMat $normalMethDir $genelistAnno $outMethNormal
# fi
# if [ ! -f ${outMethTumor}_diffMeth.mat ]; then 
#   $RSCRIPT $srcDir/model/step1-5_getDiffMethy.r --tumor $outMethTumor --normal $outMethNormal 
# fi


# geneSamplelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt
# if [ ! -f $geneSamplelist ]; then
#   $RSCRIPT $srcDir/model/step1-6_getCNVMethFreeSample.r --cnv $cnvmat --meth ${outMethTumor}_diffMeth.mat --out $geneSamplelist 
# fi

#--- DEG
tumorExp=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/brca_exp_level3_02042014.mat
normalExp=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/brca_expNormal_level3_02042014.mat
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list.geneSingleStartEnd
genesample=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt
outgsdeg=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt
# $RSCRIPT $srcDir/model/step1-7_getDEG.r --tumor $tumorExp --normal $normalExp --gene $genelist --genesample $genesample --out outgsdeg 


#-----------regulator
cernetBrca=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
out=${outgsdeg}_regulator.txt
####----stop here for develop!!!
# $PYTHON $srcDir/model/step1-8_extractCeRNETRegulator.py -i $outgsdeg -d $cernetBrca -o $out

##-----------# #----------##-----------#2nd run for nature list done 
##-----------# #----------##-----------# #-----------se----------------------------------c
##-----------# #----------##-----------# #-----------se-----------------------------c
#----------##-----------sec -get data for model
##-----------test for targetGene ESR1
#---step1: create temp dir and extract all sample,regulator,exp,cnv,snp,som data
# $srcDir/model/getGeneData.sh ESR1
#--step1.1---summary report for ESR1 

#---step2 :fireoff regression for each regulator: different type based on mutations numbers 

# $RSCRIPT
echo $END
