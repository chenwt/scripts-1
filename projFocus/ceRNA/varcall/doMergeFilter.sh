#!/bin/bash
#By: J.He
#$-cwd
#TODO: 


srcDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/
rawVarDIR=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/
filteredDIR=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/

cd $filteredDIR
if [ ! -f $1 ]; then ln -s $rawVarDIR/$1 . ; fi

# echo "$srcDIR/step-5_doAnnovarAll.sh $1 " |qsub -l mem=8g,time=8:: -N ann_test -cwd -e .log/ -o ./log 
# $srcDIR/step-5_doAnnovarAll.sh $1  
# if [ $? == 0 ]; then
# rm $1.annovar.summary.invalid_input
# rm $1.annovar.summary.log 
# rm $1_summary.csv
# else
  # echo 'ERR: annovar annotation'
# fi

annVCF=$1.annovar.summary.genome_summary.csv.vcf
echo $annVCF
$srcDIR/step-5_filterGatkAnnovar.sh $annVCF
