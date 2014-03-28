#!/bin/bash
#! -cwd
#By: J.He
#Desp.: all genes, and samples for Mar2014 run of ceRNA project
#TODO:

gslistDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 

### after runnning all step1, get all matrix ready, do this script
tglist=$gslistDir/CG_combine_02172014.list

# cnvMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix
# tgcnv=$gslistDir/tgcnv.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $cnvMat -o $tgcnv -t uniq


# expMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-23-2014.voomNormed.matrix
# tgexp=$gslistDir/tgexp.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $expMat -o $tgexp -t uniq 

# somMatp=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix
# somMat3=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix
# somMat5=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.utr5p.Mar-20-2014.matrix

tgsomp=$gslistDir/tgsomp.temp
tgsom3=$gslistDir/tgsom3.temp
tgsom5=$gslistDir/tgsom5.temp
tgsom=$gslistDir/tgsom.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $somMatp -o $tgsomp
# # $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $somMat3 -o $tgsom3
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $somMat5 -o $tgsom5
#### awk 'NR>1&&FNR!=1' $gslistDir/tgsom*.temp |~/bin/sortxh > $tgsom
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/combineSomMatrix.py -p $tgsomp -t $tgsom3 -f $tgsom5 -o ${tgsom}_${CDT}.uniq

methdiffMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth/brca_meth_l3_tumor_Mar-24-2014.matrix_diffMeth.matrix
tgmeth=$gslistDir/tgmeth.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $tglist -d $methdiffMat -o $tgmeth -t uniq

barcodelist=$gslistDir/barcode_tumor_Mar-21-2014.list
genelist=$gslistDir/CG_target_Mar-23-2014.list
tgcnv=$gslistDir/tgcnv.temp
tgmeth=$gslistDir/tgmeth.temp
tgsom=$gslistDir/tgsom.temp_Mar-24-2014.uniq
out=$gslistDir/gslist_${CDT}
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getGslistGint.py -b $barcodelist -g $genelist -c $tgcnv -m $tgmeth -s $tgsom -o $out 
# awk -F"\t|;" 'NF>10{print $0 }' $gslistDir/gslist_Mar-24-2014_CnvMethSomFree >  $gslistDir/gslist_Mar-24-2014_CnvMethSomFree.10smapMore

tumor=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix
normal=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix

###-----------CNV meth free
gslist=$gslistDir/gslist_Mar-24-2014_CnvMethFree.10more
awk -F"\t|;" 'NF>10{print $0 }' $gslistDir/gslist_Mar-24-2014_CnvMethFree >  $gslist
out=$gslist.${CDT}
$RSCRIPT /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getDEgeneSample.r --tumor $tumor --normal $normal --gslist $gslist --out $out

gslist=$gslistDir/gslist_Mar-24-2014_CnvMethFree.10smapMore.deg_20140325.txt.10more
awk -F"\t|;" 'NF>10{print $0 }' $gslistDir/gslist_Mar-24-2014_CnvMethFree.10more.deg_20140325.txt > $gslist 
cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genStat4Gslist.py -i $gslist -c $cernet -o ${gslist}_stat

## RUN /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plot/plotStatsGslist.r on MacOS(using latex)
## report dir output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslit_gint_withSom_Mar-25-2014_stats_"


###------------- Som Meth Free 
gslist=$gslistDir/gslist_Mar-24-2014_CnvMethSomFree.10smapMore
out=${gslist}.${CDT}
$RSCRIPT /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getDEgeneSample.r --tumor $tumor --normal $normal --gslist $gslist --out $out 
awk -F"\t|;" 'NF>10{print $0 }' $gslistDir/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt >  $gslistDir/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more

gslist=$gslistDir/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more
cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt

$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genStat4Gslist.py -i $gslist -c $cernet -o ${gslist}_stat

##---viz
## RUN /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plot/plotStatsGslist.r on MacOS(using latex)
## report dir output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslit_gint_Mar-25-2014_stats_"


