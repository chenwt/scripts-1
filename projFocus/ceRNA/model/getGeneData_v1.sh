#!/bin/bash
#! -cwd
#By: J.He
#Desp: generate test data for grplasso regression model, this scripted is called by step2-2_grpLassoSNP.r 
#input: <Cancer Target gene Name>  
#output: <cnt of each file> < a temp dir with exp, cnv, snp, som, sample, regulator information data

#-----------funcStart
getGeneData(){
  gene=$1
  file=$2
  output=$3
  head -1 $file > $output
  grep -w $gene $file >> $output
}
cntRecord(){
    awk 'END{print NR-1}' $1  
}
#-----------funcEdn
# gene="ESR1"
host=`hostname`
if [[ $host == "c2b2acml10.c2b2.columbia.edu" ]]; then 
  rootd="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus"
  wd=/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model
else
  rootd="/ifs/data/c2b2/ac_lab/jh3283/projFocus"
  wd=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model
fi

gene=$1
genesample=$rootd"/result/02022014/geneSamples/brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt.deg_2014-02-24.txt_regulatorSamples"
cernet=$rootd"/other/brca_ceRNA_network.txt"
##data for target
# expdata=$rootd"/result/02022014/expression/brca_exp_level3_02042014.mat_voomed_2014-02-24.mat" ##voomed transformed after 
# cnvdata=$rootd"/result/02022014/cnv/brca_cnvTumor_level2_combinedCG_02242014.mat"
# snpdata=$rootd"/result/02022014/snp/brca_snp_KWtest_combinedCG_DEG.mat.1e-06_2014-02-25.mat"

expdata=$rootd"/result/02022014/expression/brca_exp_level3_02042014.mat_voomed_2014-02-26.mat"
cnvdata=$rootd"/result/02022014/cnv/brca_cnvTumor_level2_combinedCGDEG_Regulator_02282014_uniq.mat"
snpdata=$rootd"/result/02022014/snp/brca_snp_KWtest_combinedCGDEG_Regulator.mat.1e-06_2014-02-28.mat"
# somdata=$rootd""

cd $wd
if [[ ! -d temp-$gene ]]; then mkdir temp-$gene ; fi
if [[ ! -d temp-$gene/log ]]; then mkdir temp-$gene/log ; fi

cd $wd/temp-$gene
#-----------getData 
getGeneData $gene $cernet regulator
getGeneData $gene $genesample sample.temp 
awk -F"\t" -v g=$gene '$2==g{print $0}' sample.temp > sample
rm sample.temp
getGeneData $gene $expdata ${gene}_exp 

#-----------output count info
out="$wd/temp-$gene/stat.txt"
echo $out
cntReg=`cntRecord regulator`
cntSample=`awk -F"\t" 'NR==2{print split($3,a,";")}' sample `
echo -e "TagrgedCancerGene\t$gene\tsampleNumber\t$cntSample" > $out
echo -e "RegulatorGene\tExpression\tCNV\tSNP\tSOM" >> $out

#-----------get data for all regulators
while read line
do
  geneTemp=`echo $line |awk '{print $1}'`
  # echo $geneTemp
  getGeneData $geneTemp $expdata ${geneTemp}_exp 
  getGeneData $geneTemp $cnvdata ${geneTemp}_cnv 
  getGeneData $geneTemp $snpdata ${geneTemp}_snp 
# getGeneData $geneTemp $somdata ${geneTemp}_som
  cntExp=` awk 'END{print NR-1}' ${geneTemp}_exp`
  cntCNV=` awk 'END{print NR-1}' ${geneTemp}_cnv`
  cntSNP=` awk 'END{print NR-1}' ${geneTemp}_snp`
# cntSOM=` awk 'END{print NR-1}' ${geneTemp}_som`
  cntSOM=0
  echo -e "$geneTemp\t$cntExp\t$cntCNV\t$cntSNP\t$cntSOM" >> "$out"
done < sample

echo "SUCESS"
