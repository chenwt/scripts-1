#!/bin/bash
#! -cwd
#By: J.He
#Desp.: given all drivers from lasso model, and the mutations grouped by gene, calculated the differential expression level of each driver in mutated sample versus nonmutated sample(in genomic intact samples of this driver)

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step3_funcMutKegReg/regIntSmp
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 



reglist=$CWD/candiReg.list 
awk 'NR>1{print $1}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix > $reglist
echo "dirver list done.." 
ls -alt $reglist

cnvMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/cnv/brca_cnv_l3_tumor_Mar-23-2014.matrix.uniq.matrix
tgcnv=$CWD/tgcnv.temp
expMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-23-2014.voomNormed.matrix
tgexp=$CWD/tgexp.temp
methdiffMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth/brca_meth_l3_tumor_Mar-24-2014.matrix_diffMeth.matrix
tgmeth=$CWD/tgmeth.temp
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $reglist -d $cnvMat -o $tgcnv -t uniq
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $reglist -d $expMat -o $tgexp -t uniq 
# $PYTHON $srcDir/processData/extractTargetMatFromAllMat.py  -i $reglist -d $methdiffMat -o $tgmeth -t uniq
# echo "diffmeth mat done..."
# ls -alt $tgmeth 
# echo "expression mat done..."
# ls -alt $tgexp 
# echo "cnv mat done..."
# ls -alt $tgcnv 

# ln -s /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/barcode_tumor_Mar-21-2014.list $CWD/barcode_tumor.list

barcodelist=$CWD/barcode_tumor.list
genelist=$CWD/candiReg.list
tgcnv=$CWD/tgcnv.temp
tgmeth=$CWD/tgmeth.temp
out=$CWD/gslist_${CDT}

# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getGslistGint.py -b $barcodelist -g $genelist -c $tgcnv -m $tgmeth -o $out 
# echo "genomic intact sample done..."
# head -5 $CWD/gslist_Apr-23-2014_CnvMethFree

tumor=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix
normal=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix

###-----------CNV meth free
gslist=$CWD/gslist_Apr-23-2014_CnvMethFree
awk -F";" 'NR==1||NF>1{print $0}' $gslist  >  $gslist.nonzeor
# out=$CWD/reglist.diff
# $RSCRIPT /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getDEgeneSample.r --tumor $tumor --normal $normal --gslist $gslist --out $out

# gslist=$CWD/gslist_Mar-24-2014_CnvMethFree.10smapMore.deg_20140325.txt.10more
# awk -F"\t|;" 'NF>10{print $0 }' $CWD/gslist_Mar-24-2014_CnvMethFree.10more.deg_20140325.txt > $gslist 
# cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genStat4Gslist.py -i $gslist -c $cernet -o ${gslist}_stat

# ## RUN /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plot/plotStatsGslist.r on MacOS(using latex)
# ## report dir output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslit_gint_withSom_Mar-25-2014_stats_"


# ###------------- Som Meth Free 
# gslist=$CWD/gslist_Mar-24-2014_CnvMethSomFree.10smapMore
# out=${gslist}.${CDT}
# $RSCRIPT /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getDEgeneSample.r --tumor $tumor --normal $normal --gslist $gslist --out $out 
# awk -F"\t|;" 'NF>10{print $0 }' $CWD/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt >  $CWD/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more

# gslist=$CWD/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more
# cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt

# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genStat4Gslist.py -i $gslist -c $cernet -o ${gslist}_stat

# ##---viz
# ## RUN /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plot/plotStatsGslist.r on MacOS(using latex)
# ## report dir output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslit_gint_Mar-25-2014_stats_"


