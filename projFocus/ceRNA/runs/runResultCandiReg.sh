#!/bin/bash
#$ -cwd
#By: J.He
#TODO: 
#Desp: this is the running file for all coding testing in this folder
##run on selected know BRCA genes

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
candiRegDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg

##----pickle dump data
cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpCernet.py brca_ceRNA_network.txt 

refseqTsstse=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
expTum=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
expNorm=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpExp.py $expTum $refseqTsstse 
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpExp.py $expNorm $refseqTsstse 


qsubgKR_gene() {
    gene=$1
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
    geneAnnofile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
    out=$candiRegDir/${gene}_candidateRegs_${CDT}.txt
    if [ ! -f $out ]; then 
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
    	# $cmd 
    	echo $cmd |qsub -l mem=8g,time=6:: -N ${gene}.KeyReg -e $candiRegDir/log -o $candiRegDir/log -cwd  >> qsubGKR.logs
    	tail -1 qsubGKR.logs
    else
	echo $out" existed " 	
    fi 
}


qsubGetKeyRegs() {
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
    # if [ ! -f $out ]; then
    runFlag=`grep -w $gene pid_running.txt|awk 'END{print NR}'`
    # if [ ! -f $out ] && [ $runFlag -gt 0 ] ; then
    if [ ! -f $oldout ] && [ $runFlag -eq 0 ] ; then
	# cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	# echo $cmd |qsub -l mem=16g,time=120:: -N ${gene}.KeyReg -e $CWD/log -o $CWD/log -cwd  >> $CWD/qsubGKR.logs

	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg_v4.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
	echo $cmd |qsub -l mem=8g,time=40:: -N ${gene}.KeyReg -e $CWD/log -o $CWD/log -cwd  >> $CWD/qsubGKR.logs
	tail -1 $CWD/qsubGKR.logs
	((cnt=cnt+1)) 	
    # else
	# echo -e " output existed" 
    fi
  done < $1
  echo $cnt "submitted!"
}

localgKR_gene() {
    gene=$1
    if [ ! -d $candiRegDir/log ] ; then mkdir $candiRegDir/log ; fi 
    gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
    geneAnnofile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
    out=$candiRegDir/${gene}_candidateRegs_${CDT}.txt
    if [ ! -f $out ]; then 
	cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
    	$cmd 
    	# echo $cmd |qsub -l mem=8g,time=6:: -N ${gene}.KeyReg -e $candiRegDir/log -o $candiRegDir/log -cwd  >> qsubGKR.logs
    	# tail -1 qsubGKR.logs
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
###----

# awk 'NR>1{print $1}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list > gslist.Gintset_Mar31.txt
# head -10 gslist.Gintset_Mar31.txt > input_qsubGKR_1.txt
##---test---
# qsubgKR_gene NRIP1  
# tarExp=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/test/DNMT3A_reg_exp_tumor
# $RSCRIPT $srcDir/model/step2-2_regKeyRegulators.r $tarExp 

# qsubGetKeyRegs input_qsubGKR_1.txt
# qsubgKR_gene DLG1
# qsubgKR_gene RYK 
# ~/bin/splitByN gslist.Gintset_Mar31.txt 10 
# localGetKeyRegs gslist.Gintset_Mar31.txt_2
# qsubGetKeyRegs gslist.Gintset_Mar31.txt_2 
# qsubGetKeyRegs gslist.Gintset_Mar31.txt_3 
# qsubGetKeyRegs gslist.Gintset_Mar31.txt_1 
# for i in `seq 4 10`
# for i in `seq 11 15`
# for i in `seq 16 36`
# sleep 150m
# for i in `seq 37 47`
# do
  # qsubGetKeyRegs gslist.Gintset_Mar31.txt_${i}
# done

# sleep 150m
# for i in `seq 48 60`
####-------------------------
# qst |awk 'NR>2&&$3!="QRLOGIN"{split($3,a,".");print a[1]}' > pid_running.txt
# for i in `seq 1 60` ##resub with longer time
for i in `seq 3 60` ##resub with longer time
do
  qsubGetKeyRegs gslist.Gintset_Mar31.txt_${i}
done


