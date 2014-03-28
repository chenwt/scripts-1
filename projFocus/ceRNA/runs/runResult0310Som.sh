#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 

CWD="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/som"
CDT=`~/bin/dt`
PYTHON="/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python"
srcDir="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA"

# cernetBrca=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
# genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/CG_combine_02172014.list
# out=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_allRegulator.txt
# $PYTHON $srcDir/model/step1-8_extractCeRNETRegulator.py -i $genelist -d $cernetBrca -o $out
# outall=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/GC_allGene.txt
# cat $genelist $out > $outall
# $PYTHON $srcDir/annotGeneByStartEndPos.py -i $outall -o $outall.tss 



# vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu
# genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/GC_allGene.txt.tss
# output=$CWD/brca_somTumor_allGene_2k_${CDT}.maf
# if [ ! -f $output ]; then
    # $PYTHON ~/scripts/projFocus/ceRNA/processData/getMaf.py -d $vcfDir -g $genelist -l tss -k maf -c 2000 -o $output > getAllMaf_v3.log &
    # echo "$PYTHON ~/scripts/projFocus/ceRNA/processData/getMaf.py -d $vcfDir -g $genelist -l tss -k maf -c 2000 -o $output > getAllMaf_v3.log"|qsub -l mem=16g,time=14:: -N getMaf -o $CWD/log/ -e $CWD/log/ -cwd >> $CWD/qsub.log
    # tail -1 $CWD/qsub.log
# fi

###----geting matrix....
getBarcode(){
  #input vcfDir, barcode file
  vcfDir=$1
  output=$2
  echo -n "" > $2
  for f in `ls $vcfDir/*vcf`
  do
      awk -F"\t" '$1~/^#CHROM/{print $10;exit}' $f >>$2
  done
}

barcode=$CWD/brca_somTumorWU_barcode.list
getBarcode /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu $barcode

input=$CWD/brca_somTumor_allGene_2k_Mar-11-2014.maf
dist=0
output=${input}_${CDT}.mat

##two--step
###doing it by two step..

#---1 split
input=$CWD/brca_somTumor_allGene_2k_Mar-11-2014.maf
# genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/brca_gint_target.list
output=$input.mat
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_allRegulator.txt
output=${input}_allRegulator.mat
$PYTHON ~/bin/splitByKey.py -i $input -s 800  -k $genelist -o $output-temp  

# echo "$PYTHON ~/bin/splitByKey.py -i $input -s 800  -k $genelist -o $output-temp "|qsub -l mem=10g,time=12:: -N sGene_$cdt -e $CWD/log -o $CWD/log -cwd >> qsub.log
# tail -1 qsub.log

# #---2 get matrix
barcode=$CWD/brca_somTumorWU_barcode.list
splitDir=$output-temp
dist=0
getMutMatByRegion=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getMutMatByRegion_ing.py
# output=$CWD/brca_somTumorWU.mat_gint_2k_dist0_${CDT}.matrix
output=$CWD/brca_somTumorWU.mat_allRegulator_2k_dist0_${CDT}.matrix
~/tools/python/Python_current/python $getMutMatByRegion -i $splitDir -d $dist -s $barcode -o $output 
# echo "~/tools/python/Python_current/python $getMutMatByRegion -i $splitDir -d $dist -s $barcode -o $output "|qsub -l mem=10g,time=8:: -N get2kMut -o $CWD/log -e $CWD/log -cwd

