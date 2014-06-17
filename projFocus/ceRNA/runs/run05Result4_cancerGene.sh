#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

CDT=`~/bin/dt`
PYTHON="/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python"
srcDir="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA"

sigtestDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest


##--- get target cancer gene
function grepCancerGene {
  CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/
  cglist=$CWD/CG_target_Mar-23-2014.list
  
  fgreedy=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary.netfile.txt
  flasso=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/keyRegSummary_allRunMay_05212014_0.01
  ~/bin/grepf2f $cglist $flasso $CWD/CG_target_Mar-23-2014.list.tarcg.lasso
  ~/bin/grepf2f $cglist $fgreedy $CWD/CG_target_Mar-23-2014.list.tarcg.greedy
}

##--- cancer gene ceRNA driver list
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/cancerGene
cd $CWD

flasso=$CWD/CG_target_Mar-23-2014.list.tarcg.lasso
fgreedy=$CWD/CG_target_Mar-23-2014.list.tarcg.greedy

awk '{print $2}' $flasso |tr ";" "\n" | sort |uniq > $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq
awk '{print $2}' $fgreedy |tr ";" "\n" | sort |uniq > $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq

##--- get cancer gene ceRNA dirver mutations in regions 
function regSpecDriver {
  flasso=$CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq
  fgreedy=$CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq
  mut2k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.promoter2k.Jun-16-2014.matrix
  mut1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.promoter1k.Jun-16-2014.matrix
  utr3p=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.utr3p.Jun-16-2014.matrix
  allmut=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix  

  foutlassoall=${flasso}_Jun-16-2014.all
  foutlasso2k=${flasso}_Jun-16-2014.2k
  foutlasso1k=${flasso}_Jun-16-2014.1k
  foutlasso3putr=${flasso}_Jun-16-2014.utr3p
  
  foutgreedyall=${fgreedy}_Jun-16-2014.all
  foutgreedy2k=${fgreedy}_Jun-16-2014.2k
  foutgreedy1k=${fgreedy}_Jun-16-2014.1k
  foutgreedy3putr=${fgreedy}_Jun-16-2014.utr3p
  
  ~/bin/grepf2f $flasso $mut2k $foutlasso2k
  ~/bin/grepf2f $flasso $mut1k $foutlasso1k
  ~/bin/grepf2f $flasso $utr3p $foutlasso3putr
  ~/bin/grepf2f $flasso $allmut $foutlassoall 
  
  ~/bin/grepf2f $fgreedy $mut2k $foutgreedy2k
  ~/bin/grepf2f $fgreedy $mut1k $foutgreedy1k
  ~/bin/grepf2f $fgreedy $utr3p $foutgreedy3putr
  ~/bin/grepf2f $fgreedy $allmut $foutgreedyall 

}
regSpecDriver

##---- counting numbers 
function getMutCounts {
  echo "======point mutation level" 
  echo -e "2k\t1k\tutr3\t2kutr\t1kutr\ttlasso\ttgreedy\n"
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
 
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  
  wc -l $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.all 
  wc -l $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all 
  echo "======gene level"
  echo -e "2k\t1k\tutr3\t2kutr\t1kutr\ttlasso\ttgreedy\n"

  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 5 |sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 5 |sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5 |sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5 |sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5 |sort |uniq |wc -l
 
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 5|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 5|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5|sort |uniq |wc -l
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5|sort |uniq |wc -l

  cat $CWD/cg_ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.all |cut -f 1|sort|uniq|wc -l 
  cat $CWD/cg_ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all |cut -f 1|sort|uniq|wc -l 
  
}

getMutCounts

