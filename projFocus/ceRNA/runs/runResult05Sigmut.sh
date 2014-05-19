#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
sigMutDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/
#---step1-- genrate target-keyregulator file
candiRegDataDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/data/
output=$sigMutDir/keyRegSummary_donejob_05172014
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-3_genTar_drCeRNAReglist.py -i $candiRegDataDir -o $output

# tgeneKeyRegFile=$sigMutDir/keyRegSumfile.cancergene.05052014
# file1=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summaryCG_May5/target_ceRNADriver_0.01.genelist
# file2=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summary/tgeneKeyreg_0.01
# # awk -F'\t' 'NR==FNR{a[$1]=$0;next}{print a[$1]}'  $file2 $file1 |sort|uniq> $tgeneKeyRegFile
# # awk -F'"\t' 'NR==FNR{a[$1]=$0;next}{print a[$1]}'  $file2 $file1> $tgeneKeyRegFile

expfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
mutfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero
keyRegfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/keyRegSummary_donejob_05172014_0.01
intactSmplist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list
mycode=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.py

t=fix
nrand=1000
s=max
outdir=$sigMutDir/${t}_${s}_${nrand}
if [[ ! -d $outdir ]]; then 
  mkdir $outdir
  mkdir $outdir/log
fi

output=$outdir/optCorr.result
logdir=$outdir/log/

cmd="$PYTHON $mycode -i $intactSmplist -k $keyRegfile -m $mutfile -e $expfile -t $t -r $nrand -s $s -l $logdir -o $output > $PWD/qsublog_${CDT}_${t}_${nrand}_${s}"
# echo $cmd
$cmd

##--test--
# echo "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r --vanilla --exp /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp --mut /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp --output test.2.qsub.out --ttol fix --tsel max"| qsub -l mem=4g,time=2::  -e ./log -o ./log -cwd -N test2bash

# qsub -l mem=4g,time=2:: -S /ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript -e ./log -o ./log -cwd -N test1qsub /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r --vanilla --exp /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp --mut /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp --output test.2.qsub.out --ttol fix --tsel max

#python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.py -i /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list -k /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/test/keyRegSummary.test -m /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero -e /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix -t fix -r 1000 -s max -l /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/test/log/ -o $PWD/optCorr.result

###check results in test dir
# find . -name "test_optCorr_05062014*tsv" -type f -mtime -1 -exec cat {} \;|sort -k 2|less
