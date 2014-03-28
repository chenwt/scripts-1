#!/bin/bash
#! -cwd
#By: J.He
#Desp.: identify data for all regulators
#TODO:

PYTHON=~/tools/python/Python_current/python

##-----meth
preInputMatfile(){
    ls *BRCA* >> brca_meth_level3_file.temp
    cnt27=`ls *BRCA.HumanMethylation27*sdrf.txt 2>/dev/null |wc -l`
    cnt450=`ls *BRCA.HumanMethylation450*sdrf.txt 2>/dev/null |wc -l`
    echo -n "" > input_temp
    if [ $cnt27 -gt 0 ]; then
      awk -F'\t' 'NR>1{print $28"\t"$27}' *BRCA.HumanMethylation27*sdrf.txt >> input_temp
    fi
    if [ $cnt450 -gt 0 ]; then
      awk -F'\t' 'NR>1{print $28"\t"$27}' *BRCA.HumanMethylation450*sdrf.txt >> input_temp
    fi
    grep -f brca_meth_level3_file.temp input_temp > input_softlink.temp 
}

genMethMat(){
   ##changes at 02172014
   methDir=$1
   genelist=$2
   output=$3
   cd $methDir
   preInputMatfile
   wc -l input_softlink.temp
   ~/scripts/projFocus/ceRNA/step1-2.1_softlinkFiles.sh input_softlink.temp
   $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-5_makeMethMatlevel3.py -f step1-2.1_softlinkFiles.sh.log  -g $genelist -o $output.temp 
   ~/bin/trfile $output.temp $output
   rm *temp
   rm $output.temp
   rmlns $methDir
}

methDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth
cdt=`date|awk '{print $2"-"$3"-"$6}'`
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulator.txt.tss
awk 'BEGIN{OFS=FS="\t"}
    {print $1,$2,$3-1000000,$4+1000000,$5}' $genelist > $methDir/temp

genelist=$methDir/temp
outputTumor=$methDir/brca_methTumor_Gint_Regulators_${cdt}.mat
# methDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/methlevel3/
# echo "genMethMat $methDir $genelist $outputTumor" | qsub -N getMethMat_t -l mem=10g,time=20:: -e $methDir/log -o $methDir/log -cwd

outputNormal=$methDir/brca_methNormal_Gint_Regulators_${cdt}.mat
methDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/methlevel3Normal/
# genMethMat $methDir $genelist $outputNormal

echo "genMethMat $methDir $genelist $outputNormal" | qsub -N getMethMat_n -l mem=10g,time=20:: -e $methDir/log -o $methDir/log -cwd

# methDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth
# rm $methDir/temp

##som

