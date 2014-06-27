#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

##parepare input files for ~/scripts/projFocus/ceRNA/funcfilter/step5-1_mutInBS_intSmp.py

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/cg/intSmpCnt

#function prepareInput {
  #---1
  awk 'NR==1{printf "point\ttssgene1\tchr\ttssps\ttsspe\tmutg\tmutps\tmutpe\t";for(i=5;i<=NF;i++) printf $i"\t"}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som/brca_tcga_som.utr3p.matrix.sorted.uniq > $CWD/header_mut_matrix.allsamples 
  
  fmutmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/cg/cg_ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06172014
  cat $CWD/header_mut_matrix.allsamples $fmutmat|awk '{FS=OFS="\t"}{printf $6"\t"$3;for(i=7;i<=NF;i++){printf "\t"$i};printf "\n"} '|~/bin/sortxh |uniq > $CWD/cg_ceRNA_driver_greedy_utr3p.mirBindSite.06242014.mutmat 
  
  awk 'NR==1{printf "point\ttssgene1\tchr\ttssps\ttsspe\tmutg\tmutps\tmutpe\t";for(i=5;i<=NF;i++) printf $i"\t"}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som/brca_tcga_som.promoter2k.sorted.matrix.sorted.uniq > $CWD/header_mut_matrix.allsamples 
  
  fmutmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/cg/cg_ceRNA_driver_greedy_1k.TFBindSite.mut.point_06172014


  cat $CWD/header_mut_matrix.allsamples $fmutmat |awk '{FS=OFS="\t"}{printf $6"\t"$3;for(i=7;i<=NF;i++){printf "\t"$i};printf "\n"} ' |~/bin/sortxh|uniq > $CWD/cg_ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
#}

PYTHON=~/tools/python/Python_current/python
scriptGetMutCnt=~/scripts/projFocus/ceRNA/funcfilter/step5-1_mutInBS_intSmp.py

scriptGetMutCnt=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/funcfilter/step5-1_mutInBS_intSmp_geneLevel.py

#function getCnt { 
  fmutmatrix=$CWD/cg_ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
  ftarIntSmple=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list
  ftarDriReg=$CWD/optCorr.result_flex_max_1000.tsv.significant.tarDriReg
  
  # awk 'BEGIN{FS=OFS="\t"}{print $1,$12}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary > $ftarDriReg 
  # fout=${fmutmatrix}.MutationSmpCnt
  fout=${fmutmatrix}.MutgeneSmpCnt
  # $PYTHON $scriptGetMutCnt -m ${fmutmatrix} -g $CWD/test.tarDriReg -s ${ftarIntSmple} -o ${fout} 
  $PYTHON $scriptGetMutCnt -m ${fmutmatrix} -g ${ftarDriReg} -s ${ftarIntSmple} -o ${fout} 
  ~/bin/sortxh $fout|uniq > $fout.uniq
  
  fmutmatrix=$CWD/cg_ceRNA_driver_greedy_utr3p.mirBindSite.06242014.mutmat
  # fout=${fmutmatrix}.MutationSmpCnt
  fout=${fmutmatrix}.MutgeneSmpCnt
  $PYTHON $scriptGetMutCnt -m ${fmutmatrix} -g ${ftarDriReg} -s ${ftarIntSmple} -o ${fout} 
  ~/bin/sortxh $fout |uniq > $fout.uniq
#} 

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/allg/intSmpCnt
function prepareInput {
  #---1
  awk 'NR==1{printf "point\ttssgene1\tchr\ttssps\ttsspe\tmutg\tmutps\tmutpe\t";for(i=5;i<=NF;i++) printf $i"\t"}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som/brca_tcga_som.utr3p.matrix.sorted.uniq > $CWD/header_mut_matrix.allsamples 
  
  fmutmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/allg/ceRNA_driver_greedy_utr3p.mirBindSite.mut.point_06232014
  cat $CWD/header_mut_matrix.allsamples $fmutmat |awk '{FS=OFS="\t"}{printf $6"\t"$3;for(i=7;i<=NF;i++){printf "\t"$i};printf "\n"} '|~/bin/sortxh |uniq > $CWD/ceRNA_driver_greedy_utr3p.mirBindSite.06232014.mutmat 
  
  awk 'NR==1{printf "point\ttssgene1\tchr\ttssps\ttsspe\tmutg\tmutps\tmutpe\t";for(i=5;i<=NF;i++) printf $i"\t"}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som/brca_tcga_som.promoter2k.sorted.matrix.sorted.uniq > $CWD/header_mut_matrix.allsamples 
  
  fmutmat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/allg/ceRNA_driver_greedy_1k.TFBindSite.mut.point_06232014
  cat $CWD/header_mut_matrix.allsamples $fmutmat |awk '{FS=OFS="\t"}{printf $6"\t"$3;for(i=7;i<=NF;i++){printf "\t"$i};printf "\n"} '|~/bin/sortxh |uniq> $CWD/ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
}

PYTHON=~/tools/python/Python_current/python
# scriptGetMutCnt=~/scripts/projFocus/ceRNA/funcfilter/step5-1_mutInBS_intSmp.py

scriptGetMutCnt=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/funcfilter/step5-1_mutInBS_intSmp_geneLevel.py

function getCnt { 
  fmutmatrix=$CWD/ceRNA_driver_greedy_1k.TFBindSite.06242014.mutmat
  ftarIntSmple=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list
  ftarDriReg=$CWD/optCorr.result_flex_max_1000.tsv.significant.tarDriReg
  
  # awk 'BEGIN{FS=OFS="\t"}{print $1,$12}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv.significant.summary > $ftarDriReg 

  # fout=${fmutmatrix}.MutationSmpCnt
  # $PYTHON $scriptGetMutCnt -m ${fmutmatrix} -g $CWD/test.tarDriReg -s ${ftarIntSmple} -o ${fout} 

  fout=${fmutmatrix}.MutgeneSmpCnt
  $PYTHON $scriptGetMutCnt -m ${fmutmatrix} -g ${ftarDriReg} -s ${ftarIntSmple} -o ${fout} 
  ~/bin/sortxh $fout|uniq > $fout.uniq
  
  fmutmatrix=$CWD/ceRNA_driver_greedy_utr3p.mirBindSite.06232014.mutmat
  # fout=${fmutmatrix}.MutationSmpCnt
  fout=${fmutmatrix}.MutgenenSmpCnt
  $PYTHON $scriptGetMutCnt -m ${fmutmatrix} -g ${ftarDriReg} -s ${ftarIntSmple} -o ${fout} 
  ~/bin/sortxh $fout |uniq > $fout.uniq
} 

