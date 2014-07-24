#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
###--functions 
function test_one_gene {
    t=fix
    nrand=1000
    s=max
    nperm=100
    outdir=$CWD/test
    output=$outdir/optCorr.result
    logdir=$outdir

    awk 'NR==2' $keyRegfile  > $keyRegfile.test
    keyRegfile=${keyRegfile}.test

    cmd="$PYTHON $mycode -i $intactSmplist -k $keyRegfile -m $mutfile -e $expfile -t $t -r $nrand -s $s -l $logdir -o $output > $PWD/qsublog_${CDT}_${t}_${nrand}_${s}"
    echo $cmd
    $cmd 
}


function qsub_flex_max {
    keyRegfile=$1
    # genefile=$1
    # ~/bin/grepf2f $genefile $keyRegfile $keyRegfile.test
    # keyRegfile=${keyRegfile}.test

    t=flex
    nrand=1000
    s=max
    nperm=$2

    outdir=$CWD/${t}_${s}_${nrand}

    if [[ ! -d $outdir ]]; then 
      mkdir $outdir
      mkdir $outdir/log
    fi

    output=$outdir/optCorr.result
    logdir=$outdir/log/
    cmd="$PYTHON $mycode -i $intactSmplist -k $keyRegfile -m $mutfile -e $expfile -t $t -r $nrand -p $nperm -s $s -l $logdir -o $output > $PWD/qsublog_${CDT}_${t}_${nrand}_${s}"
    echo $cmd
    # $cmd 

    echo "[START time]"`date` >> $CWD/logs.run
    $cmd >> $CWD/logs.run
    tail -1 $CWD/logs.run
    echo "[END time]"`date` >> $CWD/logs.run
}

function qsub_fix_max {
    keyRegfile=$1
    t=fix
    nrand=1000
    s=max
    nperm=100

    outdir=$CWD/${t}_${s}_${nrand}

    if [[ ! -d $outdir ]]; then 
      mkdir $outdir
      mkdir $outdir/log
    fi

    output=$outdir/optCorr.result
    logdir=$outdir/log/
    cmd="$PYTHON $mycode -i $intactSmplist -k $keyRegfile -m $mutfile -e $expfile -t $t -r $nrand -p $nperm -s $s -l $logdir -o $output > $PWD/qsublog_${CDT}_${t}_${nrand}_${s}"
    # echo $cmd
    echo "[START time]"`date` >> $CWD/logs.run
    $cmd >> $CWD/logs.run
    tail -1 $CWD/logs.run
    echo "[END time]"`date` >> $CWD/logs.run
}
##--function--end


##--function--end

#---step1--Initi

###change to tf file
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/sigReg/runJuly17
intactSmplist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist/brca_tcga_rnaseq851_aracne.tflistallsample.gslist
keyRegfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/candiReg/runJuly15/summary/tf_aracne_keyReg_runJuly15_0.01.hasExp

expfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
mutfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero

mycode=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr_permuAll.py


##---step1.2--split input file
# ~/bin/splitByN $keyRegfile 10

##--step 2--qsub-jobs

# test_one_gene
# i=0
# qsub_flex_max ${keyRegfile}_${i} 1000 &
# i=1
# qsub_flex_max ${keyRegfile}_${i} 100


# for i in `seq 3 20`
# do
#   qsub_flex_max ${keyRegfile}_${i} 100
#   sleep 180m
# done 

# i=2
# qsub_flex_max ${keyRegfile}_${i} 100

# i=3
# qsub_flex_max ${keyRegfile}_${i} 100

# i=4
# qsub_flex_max ${keyRegfile}_${i} 100

# i=1stg
# qsub_fix_max ${keyRegfile}_${i} 100

# for i in `seq 0 3`
# do
#   qsub_fix_max ${keyRegfile}_${i} 100
#   sleep 180m
# done 


# for genefile in `ls $CWD/test.gene_*`
# do 
#   echo $genefile 
#   # qsub_flex_max $genefile  
# done 


##--------------clean results

function getSumTSV {
  cd $CWD/flex_max_1000
  if [ ! -d permSummary ]; then
    mkdir permSummary
  fi 
  keyRegfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/candiReg/runJuly15/summary/tf_aracne_keyReg_runJuly15_0.01.hasExp
  for inputf in ${keyRegfile}_1stg ${keyRegfile}_0 ${keyRegfile}_1 ${keyRegfile}_2 ${keyRegfile}_3
  do 
    for g in `cut -f1 $inputf `
    do
      # awk 'NR==1||FNR!=1' optCorr.result_${g}_flex_max_1000* > permSummary/summary_${g}_result_flex_max_1000.tsv
      # awk 'FNR!=1' optCorr.result_${g}_*_flex_max_1000.permuAll_* >> permSummary/summary_${g}_result_flex_max_1000.tsv

      awk 'NR==1||FNR!=1' optCorr.result_${g}_fix_max_1000* > permSummary/summary_${g}_result_fix_max_1000.tsv
      awk 'FNR!=1' optCorr.result_${g}_*_fix_max_1000.permuAll_* >> permSummary/summary_${g}_result_fix_max_1000.tsv
    done
  done
}

# getSumTSV $CWD/flex_max_1000


# sumResult

##------get significant ones and generate plot
# run /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptSummary.r 


function singleTest {
  ##--test--
  # echo "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r --vanilla --exp /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp --mut /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp --output test.2.qsub.out --ttol fix --tsel max"| qsub -l mem=4g,time=2::  -e ./log -o ./log -cwd -N test2bash
  
  # qsub -l mem=4g,time=2:: -S /ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript -e ./log -o ./log -cwd -N test1qsub /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.r --vanilla --exp /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp --mut /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp --output test.2.qsub.out --ttol fix --tsel max
  
  #python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr.py -i /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list -k /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/test/keyRegSummary.test -m /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero -e /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix -t fix -r 1000 -s max -l /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/test/log/ -o $PWD/optCorr.result
  
  # /ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_greedyOptCorr_underwork.r --vanilla --exp /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp --mut /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp --output test.modifiedFlex --ttol flex --tsel max --nrand 10
  
  ###check results in test dir
  find . -name "test_optCorr_05062014*tsv" -type f -mtime -1 -exec cat {} \;|sort -k 2|less
}



