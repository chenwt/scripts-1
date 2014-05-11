#!/bin/bash
#$ -cwd
#By: J.He
#Desp: this is the running file for all coding testing in this folder
##run on selected know BRCA genes

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA

##----pickle dump data
cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt

refseqTsstse=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
expTum=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
expNorm=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpExp.py $expTum $refseqTsstse 
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpExp.py $expNorm $refseqTsstse 
gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list
geneAnnofile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv

qsubgKR_gene() {
    gene=$1
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    out=$candiRegDir/${gene}_candidateRegs_${CDT}.txt
    nperm=1000
    if [ ! -f $out ]; then 
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py -c $cernet -p $nperm -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
    	echo $cmd |qsub -l mem=8g,time=6:: -N ${gene}.KeyReg -e $candiRegDir/log -o $candiRegDir/log -cwd  >> $candiRegDir/qsubGKR.logs
    	tail -1 $candiRegDir/qsubGKR.logs
    else
	echo $out" existed " 	
    fi 
}


qsubGetKeyRegsSmall() {
  dirName=`echo $1|awk -F "/" '{print $NF}'`
  CWD=$candiRegDir/temp-$dirName
  echo $CWD
  if [ ! -d $CWD ] ; then mkdir $CWD ; fi 
  if [ ! -d $CWD/log ] ; then mkdir $CWD/log ; fi 
  cnt=0
  while read gene
  do 
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    runFlag=0
    if [ -f $candiRegDir/pid_running.txt ] ; then
      runFlag=`grep -w $gene $candiRegDir/pid_running.txt|awk 'END{print NR}'`
    fi 

    gene=$gene
    out=$CWD/${gene}_candidateRegs.txt
    nperm=1000

    if [ ! -f $out ] && [ $runFlag -eq 0 ] ; then
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py -c $cernet -p $nperm -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	echo $cmd |qsub -l mem=8g,time=160:: -N ${gene}.KeyReg -e $CWD/log -o $CWD/log -cwd  >> $CWD/qsubGKR.logs
	tail -1 $CWD/qsubGKR.logs
	((cnt=cnt+1))	
    fi
  done < $1
  echo $cnt "submitted!"
}

qsubGetKeyRegsBig() {
  dirName=`echo $1|awk -F "/" '{print $NF}'`
  CWD=$candiRegDir/temp-$dirName
  echo $CWD
  if [ ! -d $CWD ] ; then mkdir $CWD ; fi 
  if [ ! -d $CWD/log ] ; then mkdir $CWD/log ; fi 
  cnt=0
  while read gene
  do 
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    runFlag=0
    if [ -f $candiRegDir/pid_running.txt ] ; then
      runFlag=`grep -w $gene $candiRegDir/pid_running.txt|awk 'END{print NR}'`
    fi 

    gene=$gene
    out=$CWD/${gene}_candidateRegs.txt
    nperm=100

    if [ ! -f $out ] && [ $runFlag -eq 0 ] ; then
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py -c $cernet -p $nperm -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	echo $cmd |qsub -l mem=8g,time=160:: -N ${gene}.KeyReg -e $CWD/log -o $CWD/log -cwd  >> $CWD/qsubGKR.logs
	tail -1 $CWD/qsubGKR.logs
	((cnt=cnt+1))	
    fi
  done < $1
  echo $cnt "submitted!"
}


localgKR_gene() {
    gene=$1
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    out=$candiRegDir/${gene}_candidateRegs_${CDT}.txt
    if [ ! -f $out ]; then 
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py -p 100 -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
    	$cmd 
    else
	echo $out" existed " 	
    fi 
}


localGetKeyRegs() {
  CWD=$candiRegDir/temp-$1
  if [ ! -d $CWD ] ; then mkdir $CWD ; fi 
  if [ ! -d $CWD/log ] ; then mkdir $CWD/log ; fi 
  while read gene
  do
    gene=$gene
    echo "#-----------$gene-------------------"
    gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
    expTum=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
    expNorm=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
    geneAnnofile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
    cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
    out=$CWD/${gene}_candidateRegs_${CDT}.txt
    if [ ! -f $out ]; then
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	$cmd > $CWD/${gene}.local_stdout
    else
	echo -e $out" existed\n please remove at first to redo!" 
    fi
  done < $1
}

getErrGene(){
  ls -alt temp-gslist.Gintset_Mar31.txt_*/log/*.e*|awk '$5!="514"&&$5!="582"{print $}'
}

###----------

candiRegDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/
##step-1 devide target into small and large
gslistStat=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt_stat.GintRegCount
tgSmall=$candiRegDir/tgene_small.txt
tgBig=$candiRegDir/tgene_big.txt
# awk '$2>200{print $1}' $gslistStat > $tgBig 
# awk '$2<=200{print $1}' $gslistStat > $tgSmall 

gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list
# ~/bin/splitByN $tgSmall 150 
# ~/bin/splitByN $tgBig 150 

##step-2 run grplasso to both group
## test one gene
# qsubgKR_gene BCL9 
# qsubgKR_gene ACTA2
# qsubGetKeyRegsSmall ${tgSmall}_1 

# for i in `seq 8 20 `
# do
#   qsubGetKeyRegsSmall ${tgSmall}_${i} 
#   sleep 280m
# done

# for i in `seq 1 3 `
# do
#   qsubGetKeyRegsSmall ${tgBig}_${i} 
#   sleep 280m
# done

##step2.4 check job running
#ls -1 temp-tgene_*/*txt|awk -F"-|_|/" '{print $5}' > jobs.done
# grep -v -w -f jobs.done tgene_small.txt > jobs.small.fail
# grep -w -f jobs.small.fail /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt_stat.GintRegCount |awk '$2 ==1{print $1}' > jobs.small.fail.1reg
# grep -w -f jobs.small.fail /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt_stat.GintRegCount |awk '$2>1{print $1}' > jobs.small.fail.gt1reg

# qsubGetKeyRegsSmall jobs.small.fail.gt1reg

##step-3 calculate summary 
##-------runing the 47 remaining ones with large memory and long time  
qsubGetKeyRegs() {
  # CWD=$candiRegDir/temp-bigTarget
  CWD=$candiRegDir/temp-$1
  if [ ! -d $CWD ] ; then mkdir $CWD ; fi 
  if [ ! -d $CWD/log ] ; then mkdir $CWD/log ; fi 
  cnt=0
  while read gene
  do
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    gene=$gene
    gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
    geneAnnofile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
    out=$CWD/${gene}_candidateRegs_${CDT}.txt
    oldout=$CWD/${gene}_candidateRegs_*2014.txt
    if [ ! -f $oldout ] ; then
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v5.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	# cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	echo $cmd |qsub -l mem=20g,time=140:: -N ${gene}.KeyReg -e $CWD/log -o $CWD/log -cwd  >> $CWD/qsubGKR.logs
	tail -1 $CWD/qsubGKR.logs
	((cnt=cnt+1)) 	
    fi
  done < $1
  echo $cnt "submitted!"
}

## get all cancer gene
# grep -w -f $candiRegDir/tgene_small.txt /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_target_Mar-23-2014.list > /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summaryCG/cancergene.list  

# grep -w -f $candiRegDir/tgene_big.txt /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/CG_target_Mar-23-2014.list >> /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summaryCG/cancergene.list  


# CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg
# resDir=$CWD/run-Apr-1-2014
# # mkdir $resDir 
# cd $resDir
# # cp $CWD/temp-*/*txt . 
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/getKeyRegStats.py -d $resDir -o $resDir/kegRegs_${CDT}.summary
# sort $resDir/kegRegs_${CDT}.summary.driverRegs.list |uniq > $resDir/kegRegs_${CDT}.summary.driverRegs.list.uniq



##step4---get some summary stats

# awk 'NR==FNR{a[$1]=$2;next}{print $1,a[$1],$2}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt_stat.GintRegCount candiReg_summary_05072014_v2.txt_0.01.regCount > candiReg_summary_05072014_v2.txt_0.01.regCount_GintReg_keyReg_05112014
