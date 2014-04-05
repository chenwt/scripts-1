#!/bin/bash
#$ -cwd
#By: J.He
#Desp.: model to get driver mutations 

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
candiRegDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg
sigMutDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut

##-----------group mutations by binding sites
cupidBSfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred.3cols
gene3pUTRfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_3pUTR.tsv
mutfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix
output=$sigMutDir/brca_som_utr3p_cupidBS_${CDT}.matrix
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/mapMut2BS.py -s $cupidBSfile -a $gene3pUTRfile -m $mutfile -o $output 
##-----------selected drive mutation sets


##-----------associated with expression




