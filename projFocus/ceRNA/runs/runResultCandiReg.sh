#!/bin/bash
#$ -cwd
#By: J.He
#TODO: 
#Desp: this is the running file for all coding testing in this folder
##run on selected know BRCA genes

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
candiRegDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg

##----pickle dump data
cernet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpCernet.py brca_ceRNA_network.txt 

refseqTsstse=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
expTum=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
expNorm=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpExp.py $expTum $refseqTsstse 
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpExp.py $expNorm $refseqTsstse 

gene="PTEN"
gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
geneAnnofile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
out=$candiRegDir/${gene}_candidateRegs_${CDT}
cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out"
# echo $cmd
$cmd

# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-1_getKeyReg.py -c $cernet -g $gene  -t $expTum -n $expNorm -a $geneAnnofile -l $gslist -o $out 
