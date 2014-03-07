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
END="#---[THE END]---"
ERR="[ERROR:]"
cdt=`~/bin/dt`

#------------step1-prepareData
##-----------# #----------##-----------# #-----------sec
##-----------# #----------##-----------# 
###---second run for nature 12912 genelist
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_nature12912_novel.list
cnvmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_naturePubNovalGene_02172014.mat
# if [ ! -f $cnvmat ]; then
#   getCNVMat $genelist $cnvmat 
# fi

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

##--filter somatic mutation
# vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered
vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_UCceRNETPlusNature12929.list.geneSingleStartEnd
somDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som
cdt=`~/bin/dt`
output=$somDir/brca_somTumorWU_combinedCG_$cdt.mat

# if [ ! -f $output ]; then
#     $PYTHON $srcDir/varcall/step-6_getMAF.py -d $vcfDir -g $genelist -l gene -k mut -o $output &
#     echo "$PYTHON $srcDir/varcall/step-6_getMAF.py -d $vcfDir -g $genelist -l gene -k mut -o $output" |qsub -l mem=10g,time=12:: -N getGeneMut -o $somDir/log/ -e $somDir/log/ -cwd >> $somDir/qsub.log
#     tail -1 qsub.log
# fi

# ##---
# sommat=$somDir/brca_somTumor_combinedCG_20140301.mat
# cnvmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_naturePubNovalGene_02172014.mat
# methmat=${outMethTumor}_diffMeth.mat
# geneSamplelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_combinedCG_CNVMethSomFree_20140302.txt

# methmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth/brca_methTumor_combinedCG_03032014.mat_diffMeth.mat
# sommat=$somDir/brca_somTumorWU_combinedCG_Mar-5-2014.mat.mat
# cnvmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_combinedCG_02242014.mat
# geneSamplelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_$cdt.txt
# if [ ! -f $geneSamplelist ]; then
#     $RSCRIPT $srcDir/model/step1-6_getCNVMethFreeSample.r --cnv $cnvmat --meth $methmat --som $sommat --out $geneSamplelist 
# fi

# # ##-----------deg ceRNA targets

# tumMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/brca_exp_level3_02042014.mat
# normMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/brca_expNormal_level3_02042014.mat
# gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt
# output=$gslist.deg_$cdt.txt
# if [ ! -f $output ]; then
#     $RSCRIPT $srcDir/processData/step1-2_getDEgeneSample.r --tumor $tumMat --normal $normMat --gslist $gslist --out $output 
# fi

# #-----------regulator
cernetBrca=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
# outgsdeg=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt
outgsdeg=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt
# out=${outgsdeg}_regulator.txt
# $PYTHON $srcDir/model/step1-8_extractCeRNETRegulator.py -i $outgsdeg -d $cernetBrca -o $out
# $PYTHON $srcDir/annotGeneByStartEndPos.py -i $out -o $out.tss

##-----------------------------------------
##-----------------------------------------
##-------prepare data secion---------------


##--regulator cnv


##--regulator snp


##--regulator somatic maf 
#------test------
# vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test
# if [ -f $somDir/test.out.mat ]; then  rm $somDir/test.out.mat ; fi
# output=$somDir/test.out.mat
#-----testend-----
vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt_regulator.txt.singleStartEnd
# output=$somDir/brca_somTumorWU_combinedCG_Regulator_$cdt.mat
output=$somDir/brca_somTumorWU_combinedCG_Regulator_v3_$cdt.mat
# if [ ! -f $output ]; then
    # $PYTHON $srcDir/varcall/step-6_getMAF_v3.py -d $vcfDir -g $genelist -l tss -k maf -o $output > getRegMaf_v3.log  &
    # echo "$PYTHON $srcDir/varcall/step-6_getMAF_v3.py -d $vcfDir -g $genelist -l tss -k maf -o $output " |qsub -l mem=16g,time=14:: -N newgetRegMaf -o $somDir/log/ -e $somDir/log/ -cwd >> $somDir/qsub.log
    # $PYTHON $srcDir/varcall/step-6_getMAF.py -d $vcfDir -g $genelist -l tss -k maf -o $output  &
    # echo "$PYTHON $srcDir/varcall/step-6_getMAF.py -d $vcfDir -g $genelist -l tss -k maf -o $output " |qsub -l mem=16g,time=14:: -N getRegMaf -o $somDir/log/ -e $somDir/log/ -cwd >> $somDir/qsub.log
    # tail -1 qsub.log
# fi

# vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu_rescue/
# genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt_regulator.txt.singleStartEnd
# output=$somDir/brca_somTumorWU_combinedCG_Regulator_resuce1sample_$cdt.mat
# $PYTHON $srcDir/varcall/step-6_getMAF_v3.py -d $vcfDir -g $genelist -l tss -k maf -o $output  &
# sed -c -i 1d brca_somTumorWU_combinedCG_Regulator_resuce1sample_Mar-7-2014.mat
# cat brca_somTumorWU_combinedCG_Regulator_resuce1sample_Mar-7-2014.mat  >>$somDir/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat

###----geting matrix....
getBarcode(){
  vcfDir=$1
  output=$2
  echo -n "" > $2
  for f in `ls $vcfDir/*vcf`
  do
    awk -F"\t" '$1~/^#CHROM/{print $10;exit}' $f >>$2
  done
}

barcode=brca_somTumorWU_barcode.list
# getBarcode /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu $barcode 

input=$somDir/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulator.txt.tss
dist=100
output=${input}_"dist"${dist}_$crt.mat
cmd="$PYTHON $srcDir/processData/getMutMatByRegion.py -i $input -g $genelist -s $barcode -d 100 -o $output"
echo $cmd
$cmd & 

# ##---normal sample somatic mat
# vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/normal
# genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_UCceRNETPlusNature12929.list.geneSingleStartEnd
# output=$somDir/brca_somTumor_combinedCG_normal_20140303.mat
# # if [ ! -f $output ]; then
# #     # $PYTHON $srcDir/varcall/step-6_getMAF.py -d $vcfDir -g $genelist -l gene -k mut -o $output 
#     echo "$PYTHON $srcDir/varcall/step-6_getMAF.py -d $vcfDir -g $genelist -l gene -k mut -o $output" |qsub -l mem=10g,time=12:: -N getGeneMut -o $somDir/log/ -e $somDir/log/ -cwd >> $somDir/qsub.log
#     tail -1 qsub.log
# fi


##------------------------------------------ 
##------------------------------------------

#----------##-----------sec -get data for model
##-----------test for targetGene ESR1
#---step1: create temp dir and extract all sample,regulator,exp,cnv,snp,som data
# $srcDir/model/getGeneData.sh ESR1
#--step1.1---summary report for ESR1 

#---step2 :fireoff regression for each regulator: different type based on mutations numbers 

# $RSCRIPT
echo $END
