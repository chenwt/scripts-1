#!/bin/bash
#! -cwd
#By: J.He
#Desp.: all genes, and samples for Mar2014 run of ceRNA project
#TODO:

gslistDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 

### after runnning all step1, get all matrix ready, do this script
tglist=$gslistDir/allTaget.list

# awk 'NR>1{print $1"\n"$2}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt |sort|uniq > $tglist 

cnvMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix
tgcnv=$gslistDir/tgcnv.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $cnvMat -o $tgcnv -t uniq


expMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-23-2014.voomNormed.matrix
tgexp=$gslistDir/tgexp.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $expMat -o $tgexp -t uniq 

somMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix
tgsom=$gslistDir/tgsom.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $somMat -o $tgsom

methdiffMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth/brca_meth_l3_tumor_Mar-24-2014.matrix_diffMeth.matrix
tgmeth=$gslistDir/tgmeth.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $methdiffMat -o $tgmeth -t uniq

barcodelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/barcode_tumor_Mar-21-2014.list
genelist=$gslistDir/allTaget.list
tgcnv=$gslistDir/tgcnv.temp
tgmeth=$gslistDir/tgmeth.temp
tgsom=$gslistDir/tgsom.temp
out=$gslistDir/gslist
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getGslistGint.py -b $barcodelist -g $genelist -c $tgcnv -m $tgmeth -s $tgsom -o $out 
# awk -F"\t|;" 'NF>10{print $0 }' $gslistDir/gslist_CnvMethFree >  ${out}_CnvMethFree.10smapMore
# awk -F"\t|;" 'NF>10{print $0 }' $gslistDir/gslist_CnvMethSomFree >  ${out}_CnvMethSomFree.10smapMore

tumor=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix
normal=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix

#####-------------dysregulated CNV Meth Som Free 
gslist=$gslistDir/gslist_CnvMethSomFree.10smapMore
out=${gslist}
# $RSCRIPT /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getDEgeneSample.r --tumor $tumor --normal $normal --gslist $gslist --out $out 
out=$gslistDir/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt
# awk -F"\t|;" 'NF>10{print $0 }' $out > $out.10more 


##---summary for dysreg Gint
gslist=$gslistDir/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt
cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genStat4Gslist.py -i $gslist -c $cernet -o ${gslist}_stat

# ##---viz
# ## RUN /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plot/plotStatsGslist.r on MacOS(using latex)
# ## report dir output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslit_gint_Mar-25-2014_stats_"


