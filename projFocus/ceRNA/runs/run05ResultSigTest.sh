#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/sigTest
cd $CWD


###get all regulator 
# cat keyRegSummary_allRunMay_05212014_0.01 |awk '{print $2}'|tr ";" "\n" |sort |uniq > ceRNA_driver_lasso_${CDT}.uniq
# cat optCorr.result_flex_max_1000.tsv.significant.summary.netfile.txt | awk '{print $2}'|tr ";" "\n" |sort |uniq > ceRNA_driver_greedy_${CDT}.uniq 



##--extract utr3p, promoter2k mutation


CDT=`~/bin/dt`
PYTHON="/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python"
srcDir="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA"


# python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grepf2f.py -i ceRNA_driver_lasso_Jun-13-2014.uniq -d brca_tcga_som.promoter2k.sorted.matrix.sorted.uniq -o ceRNA_driver_lasso_Jun-13-2014.uniq.promoter2k &

# python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grepf2f.py -i ceRNA_driver_greedy_Jun-13-2014.uniq -d brca_tcga_som.utr3p.matrix.sorted.uniq -o ceRNA_driver_greedy_Jun-13-2014.uniq.utr3p  &

function regSpecDriver {
  flasso=$CWD/ceRNA_driver_lasso_Jun-13-2014.uniq
  fgreedy=$CWD/ceRNA_driver_greedy_Jun-13-2014.uniq 
  mut2k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.promoter2k.Jun-16-2014.matrix
  mut1k=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.promoter1k.Jun-16-2014.matrix
  utr3p=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.utr3p.Jun-16-2014.matrix
  allmut=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix  

  foutlasso2k=${flasso}_Jun-16-2014.2k
  foutlasso1k=${flasso}_Jun-16-2014.1k
  foutlasso3putr=${flasso}_Jun-16-2014.utr3p
  foutlassoall=${flasso}_Jun-16-2014.all
   
  
  foutgreedy2k=${fgreedy}_Jun-16-2014.2k
  foutgreedy1k=${fgreedy}_Jun-16-2014.1k
  foutgreedy3putr=${fgreedy}_Jun-16-2014.utr3p
  foutgreedyall=${fgreedy}_Jun-16-2014.all
  
  ~/bin/grepf2f $flasso $mut2k $foutlasso2k
  ~/bin/grepf2f $flasso $mut1k $foutlasso1k
  ~/bin/grepf2f $flasso $utr3p $foutlasso3putr
  ~/bin/grepf2f $flasso $allmut $foutlassoall 
  
  ~/bin/grepf2f $fgreedy $mut2k $foutgreedy2k
  ~/bin/grepf2f $fgreedy $mut1k $foutgreedy1k
  ~/bin/grepf2f $fgreedy $utr3p $foutgreedy3putr
  ~/bin/grepf2f $fgreedy $allmut $foutgreedyall 

}
#regSpecDriver

function getMutCounts {
  echo "======point mutation level" 
  echo -e "2k\t1k\tutr3\t2kutr\t1kutr\ttlasso\ttgreedy\n"
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
 
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 2,6,7|sort |uniq |wc -l
  
  wc -l $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.all 
  wc -l $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all 
  echo "======gene level"
  echo -e "2k\t1k\tutr3\t2kutr\t1kutr\ttlasso\ttgreedy\n"

  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 5 |sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 5 |sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5 |sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5 |sort |uniq |wc -l
  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5 |sort |uniq |wc -l
 
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k |cut -f 5|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k |cut -f 5|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.2k $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5|sort |uniq |wc -l
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.1k $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.utr3p |cut -f 5|sort |uniq |wc -l

  cat $CWD/ceRNA_driver_lasso_Jun-13-2014.uniq_Jun-16-2014.all |cut -f 1|sort|uniq|wc -l 
  cat $CWD/ceRNA_driver_greedy_Jun-13-2014.uniq_Jun-16-2014.all |cut -f 1|sort|uniq|wc -l 
  
}

getMutCounts
