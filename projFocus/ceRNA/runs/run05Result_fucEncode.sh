#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-2_dnaseSite

function getBedfile {
  fgallmut=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
  fcgmut=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all
  
  fgallmbed=$CWD/brca_greedy_ceRNAdirverMutation.all_06182014.bed
  fcgmbed=$CWD/cg_brca_greedy_ceRNAdirverMutation.all_06182014.bed
  
  awk 'BEGIN{FS=OFS="\t"}{print "chr"$2,$3,$4,$1}' $fgallmut >$fgallmbed 
  awk 'BEGIN{FS=OFS="\t"}{print "chr"$2,$3,$4,$1}' $fcgmut > $fcgmbed
  
}

fgallmbed=$CWD/brca_greedy_ceRNAdirverMutation.all_06182014.bed
fcgmbed=$CWD/cg_brca_greedy_ceRNAdirverMutation.all_06182014.bed
  

function getMutInDnaseMCF7 {
  dDnaseMCF7=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/encode/dnaseBS/MCF7
  cnt=0
  for fsig in `ls $dDnaseMCF7`
  do
    ((cnt=$cnt+1))
    echo "Processing $cnt "$dDnaseMCF7/$fsig
    /ifs/home/c2b2/ac_lab/rs3412/tools/bedtools/bin/bedtools intersect -wb -a $fcgmbed   -b $dDnaseMCF7/$fsig > $fcgmbed.$fsig.${CDT} &
    /ifs/home/c2b2/ac_lab/rs3412/tools/bedtools/bin/bedtools intersect -wb -a $fgallmbed -b $dDnaseMCF7/$fsig > $fgallmbed.$fsig.${CDT}  
  done
  
  ### then run /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/funcfilter/step5-2_mutInDnase_mcf7.r  to get result:q
  
}

