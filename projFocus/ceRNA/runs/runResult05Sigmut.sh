#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
sigMutDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/
#---step1-- genrate target-keyregulator file
candiRegDataDir=
output=
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-3_genTar_drCeRNAReglist.py -i $candiRegDataDir -o $output

tgeneKeyRegFile=$sigMutDir/keyRegSumfile.cancergene.05052014
file1=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summaryCG_May5/target_ceRNADriver_0.01.genelist
file2=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summary/tgeneKeyreg_0.01
# awk -F'\t' 'NR==FNR{a[$1]=$0;next}{print a[$1]}'  $file2 $file1 |sort|uniq> $tgeneKeyRegFile
# awk -F'"\t' 'NR==FNR{a[$1]=$0;next}{print a[$1]}'  $file2 $file1> $tgeneKeyRegFile

expfile=
mutfile=
out=
jobname=`echo $expfile| awk -F"." '{print $1}'`
echo " /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r --vanilla --exp $expfile --mut $mutfile --output $out"| qsub -l mem=8g,time=2:: -S /ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript -e $sigMutDir/log -o $sigMutDir/log -cwd -N  



###check results in test dir
# find . -name "test_optCorr_05062014*tsv" -type f -mtime -1 -exec cat {} \;|sort -k 2|less
