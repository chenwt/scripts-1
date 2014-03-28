#!/bin/bash
#! -cwd
#By: J.He
#Desp: generate test data for regression model to identify contributing regulators, this scripted is called by next step regression model 
##: compared to v2, extracting all CNV, SNP, Meth, and SOM data for all regulators
#input: <Cancer Target gene Name> <output dir>  <optional: ceRNET network file> <optional: expression file>, <optional: geneSample file>
#output: <cnt of each file> < a temp dir with exp, cnv, snp, som, sample, regulator information data

gene=$1
outDir=$2
#-----------funcStart
kgetGeneData(){
  gene=$1
  file=$2
  output=$3
  head -1 $file > $output
  grep -w $gene $file >> $output
}
cntRecord(){
    awk 'END{print NR-1}' $1  
}
checkStat(){
  MSG=$1
  if [[ $? != 0 ]]; then 
    echo "ERROR " $MSG
    exit
  # else
    # echo "Next"
  fi
}
#-----------funcEdn

##init
host=`hostname`
if [[ $host == "c2b2acml10.c2b2.columbia.edu" ]]; then 
  rootd="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus"
  wd=/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model
else
  rootd="/ifs/data/c2b2/ac_lab/jh3283/projFocus"
  wd=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model
fi
##setup
genesample=$rootd"/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt"
cernet=$rootd"/other/brca_ceRNA_network.txt"
expdata=$rootd"/result/02022014/expression/brca_exp_level3_02042014.mat_voomed_2014-02-26.mat"
cnvdata=$rootd"/result/02022014/cnv/brca_cnvTumor_level2_combinedCGDEG_Regulator_02282014_uniq.mat"
snpdata=$rootd"/result/02022014/snp/brca_snp_KWtest_combinedCGDEG_Regulator.mat.1e-06_2014-02-28.mat"
somdata=$rootd"/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist1000_Mar-7-2014.matrix"
# methdata=$rootd""

if [[ ! -d $outDir ]]; then 
    mkdir $outDir
fi
if [[ ! -d $outDir/log ]]; then 
  mkdir $outDir/log 
fi
checkStat "mkdir"

##---get Gint Samples
awk -F"\t|;" -v g=$gene '$1==g{for (i=2;i<=NF;i++) print $i;}' $genesample > $outDir/samples.txt
checkStat "get samples"

##---getRegulator
grep -w $gene $cernet |awk -F"\t" '{print $1"\n"$2}' |sort|uniq |awk -v g=$gene '$1!=g' > $outDir/regulators.txt
checkStat "getRegulator "

#-----------get gene expression  
head -1 $expdata > $outDir/exp.mat
awk -F"\t" -v g=$gene 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $expdata >> $outDir/exp.mat 
for line in `cat $outDir/regulators.txt`
do
  awk -F"\t" -v g=$line 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $expdata >> $outDir/exp.mat 
done
checkStat "get expression"

#--------get cnv, snp and meth
head -1 $cnvdata > $outDir/cnv.mat
awk -F"\t" -v g=$gene 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $cnvdata >> $outDir/cnv.mat 
for line in `cat $outDir/regulators.txt`
do
  awk -F"\t" -v g=$line 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $cnvdata >> $outDir/cnv.mat 
done
checkStat "get cnv"

head -1 $snpdata > $outDir/snp.mat
awk -F"\t" -v g=$gene 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $snpdata >> $outDir/snp.mat 
for line in `cat $outDir/regulators.txt`
do
  awk -F"\t" -v g=$line 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $snpdata >> $outDir/snp.mat 
done
checkStat "get snp"

###---TODO
# head -1 $methdata > $outDir/meth.mat
# awk -F"\t" -v g=$gene 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $methdata >> $outDir/meth.mat 
# for line in `cat $outDir/regulators.txt`
# do
#   awk -F"\t" -v g=$line 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $methdata >> $outDir/meth.mat 
# done
# checkStat "get meth"

#--------get somatic mutations

head -1 $somdata > $outDir/som.mat
awk -F"\t" -v g=$gene 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $somdata >> $outDir/som.mat 
for line in `cat $outDir/regulators.txt`
do
  awk -F"\t" -v g=$line 'BEGIN{FS=OFS="\t"} $1==g{print $0}' $somdata >> $outDir/som.mat 
done
checkStat "get som"


echo "SUCESS"
