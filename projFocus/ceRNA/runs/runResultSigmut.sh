#!/bin/bash
#$ -cwd
#By: J.He
#Desp.: model to get driver mutations 

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
# srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
crnsc=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
candiRegDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg
sigMutDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut

# ##-----------group mutations by binding sites
# cupidBSfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred.3cols
# gene3pUTRfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_3pUTR.tsv
# mutfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix
# output=$sigMutDir/brca_som_utr3p_cupidBS_${CDT}.matrix
# # $PYTHON $crnsc/model/mapMut2BS.py -s $cupidBSfile -a $gene3pUTRfile -m $mutfile -o $output 

##==================================================
#-----------step1 grouping 
combineMutFile=$sigMutDir/step1_grouping/brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.uniq
# $PYTHON $srcDir/model/step3-1_groupMutation.py -i $combineMutFile 
combineMutFile=$sigMutDir/step1_grouping/brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.uniq.groupByGene.binary.matrix

#-----------step2 get mutated driver

keyRegResDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/data/
sigKeyRegDir=$sigMutDir/step2_mutKeyReg/sigKeyReg
# awk -F"\t" '$5>0{print $1}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/kegRegs_Apr-10-2014.summary|sort|uniq > step2_mutKeyReg/sigCancergGene.list
keyRegSummary=$sigMutDir/step2_mutKeyReg/sigCancergGene.list

function exSigRegMutMatrix() {
  if [ ! -d $sigKeyRegDir ]; then
    mkdir $sigKeyRegDir
  fi
  
  rm -f $sigKeyRegDir/* 
  cd $sigKeyRegDir
  while read sigReg
  do 
    ln -s $keyRegResDir/${sigReg}_*.txt . 
  done < $keyRegSummary 
  cd $sigMutDir
  driverReglist=$sigMutDir/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq
  driverMutlist=$driverReglist.mut.matrix
  # $PYTHON $crnsc/model/step3-2_extractMutKeyReg.py -g $driverReglist -m $combineMutFile -o $driverMutlist 
  
}
# exSigRegMutMatrix
driverReglist=$sigMutDir/step2_mutKeyReg/keyReg_Apr-20-2014.summary.driverRegs.list
driverMutlist=$driverReglist.mut.matrix
 
# $PYTHON $crnsc/model/step3-2_extractMutKeyReg.py -g $driverReglist -m $combineMutFile -o $driverMutlist

#-----------step3 get functional mutated driver _not finalized 
##jump to step 4
function fncSigRegMut() {
  expMat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
  driverMutM=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix
  funcDriverM=$sigMutDir/step3_funcMutKegReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix.func.matrix
  $PYTHON $crnsc/model/step3-3_getFuncMutDriver.py -m $driverMutM -e $expMat -o $funcDriverM
}

#-----------step4 MSMC to get driver Set for target gene
driverMutM=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix
msmcDir=$sigMutDir/step4_MSMC
gslistfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
regmutfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step1_grouping/brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.uniq.groupByGene.binary.matrix
# out=$msmcDir/mutRegSamp.numbers4
dMutSets=$msmcDir/driverMutSet.${CDT}

#split all the sigMut into three region
sigKeyRegDir=$sigMutDir/step2_mutKeyReg/sigKeyReg
function lnsfiles(){
  smallD=$sigMutDir/step4_MSMC/small  
  bigD=$sigMutDir/step4_MSMC/big 
  if [ ! -d $smallD ] ; then  mkdir $smallD ; fi
  if [ ! -d $bigD  ] ; then mkdir $bigD ; fi 
  rm -f $smallD/*
  rm -f $bigD/*
  mins=2
  cuts=30
  cnteff=0
  while read line 
  do
    cntSmp=`echo $line | awk '{print $2}'` 
    g=`echo $line | awk '{print $1}'` 
    if [[ $cntSmp -le $cuts ]] && [[ $cntSmp -gt $mins ]] ; then
      file=`readlink -f $sigKeyRegDir/${g}_*.txt`
      if [ -f $file ]; then 
	(( cnteff=cnteff+1))
	cd $smallD
	ln -s $sigKeyRegDir/${g}_*.txt .  
      fi
    fi 
    if [[ $cntSmp -gt $cuts ]] ; then
      file=`readlink -f $sigKeyRegDir/${g}_*.txt`
      if [ -f $file ]; then
	(( cnteff=cnteff+1))
	 cd $bigD
	 ln -s $sigKeyRegDir/${g}_*.txt .  
      fi 
    fi
  done < /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/mutRegSamp.numbers.Apr172014_sig     
  echo -e "number of cancer gene has mutsample more than 2\t"$cnteff
}
# lnsfiles 

function runMSMC(){
  smallD=$sigMutDir/step4_MSMC/small  
  bigD=$sigMutDir/step4_MSMC/big 
  mutM=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix
  gsL=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
  fileD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/small
  alpha=0.85
  out=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/brca_dMutSet_msmc_small_${CDT}.tsv
  jname=msmc_`echo ${CDT}|awk -F"-" '{print $1"_"$2}'`
  # echo "$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD -o $out " |qsub -l mem=20g,time=48:: -e $sigMutDir/log -o $sigMutDir/log -N msmcV2_$jname -cwd > $sigMutDir/qsub.log
  # tail -1 $sigMutDir/qsub.log
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD  -t maxfreq -o $out & 
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD  -t minsize -o $out & 
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -m $mutM -g $gsL -d $fileD  -t minsize -o $out & 
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -m $mutM -g $gsL -d $fileD  -t maxfreq -o $out & 

  fileD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/big
  out=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/brca_dMutSet_msmc_big_${CDT}.tsv
  # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD  -t maxfreq -o $out & 
  echo "$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD -o $out -t minsize" |qsub -l mem=20g,time=48:: -e $sigMutDir/log -o $sigMutDir/log -N $jname -cwd > $sigMutDir/qsub.log
  tail -1 $sigMutDir/qsub.log
  echo "$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD -o $out -t maxfreq" |qsub -l mem=30g,time=48:: -e $sigMutDir/log -o $sigMutDir/log -N $jname -cwd > $sigMutDir/qsub.log
  tail -1 $sigMutDir/qsub.log
  # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -a $alpha -m $mutM -g $gsL -d $fileD  -t minsize -o $out &
}
runMSMC

# # echo "$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_MSMC.py -m $mutM -g $gsL -d $fileD -o $out " |qsub -l mem=20g,time=48:: -e $sigMutDir/log -o $sigMutDir/log -N msmcV1Apr18 -cwd > $sigMutDir/qsub.log
# tail -1 $sigMutDir/qsub.log


# echo $out
# echo -e "targetGene\tgintSmpCnt\tmutReg" > $out
# for keyRegfile in `ls $keyRegResDir/*txt`
# do 
#     output=$msmcDir/`echo  $keyRegfile|awk -F"/" '{split($NF,a,"_"); print a[1]}'`.msmc  
#     $PYTHON $crnsc/model/step3-4_MSMC.py -m $regmutfile -g $gslistfile -c $keyRegfile -o  $output >> $out 
# done


