#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 


#-----------function
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/runs/runPipeline_step1.sh 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/
RSCRIPT=/ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript
#-----------execute
END="#---[THE END]---"
#-step1-prepareData

genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list
##output
cnvmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_naturePubNovalGene_02172014.mat
if [ ! -f $cnvmat ]; then
  getCNVMat $genelist $cnvmat 
fi

#---
tumorMethDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/methlevel3/
normalMethDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/methlevel3Normal/
genelistAnno=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list.geneSingleStartEnd
outMethTumor=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth/brca_methTumor_nature12912NovelGene_02172014.mat
outMethNormal=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth/brca_methNormal_nature12912NovelGene_02172014.mat
if [ ! -f $outMethTumor ]; then
  genMethMat $tumorMethDir $genelistAnno $outMethTumor  
fi
if [ ! -f $outMethNormal ]; then
  genMethMat $normalMethDir $genelistAnno $outMethNormal
fi
if [ ! -f ${outMethTumor}_diffMeth.mat ]; then 
  $RSCRIPT $srcDir/model/step1-5_getDiffMethy.r --tumor $outMethTumor --normal $outMethNormal 
fi

geneSamplelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt
if [ ! -f $geneSamplelist ]; then
  $RSCRIPT $srcDir/model/step1-6_getCNVMethFreeSample.r --cnv $cnvmat --meth ${outMethTumor}_diffMeth.mat --out $geneSamplelist 
fi

#--- DEG


echo $END
