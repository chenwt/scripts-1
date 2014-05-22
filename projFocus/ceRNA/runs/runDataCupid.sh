#!/bin/bash
#$ -cwd
#By: J.He
#Desp.: 
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction
PYTHON=/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/

# awk '{print $1, $2, $4}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred | uniq > CupidPred.3cols

# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/parseBS/fasta2tsv.py refseq_hg19_refflat_cupid3pGene.fasta refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv 
# uniq refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv > refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv.uniq

seqof3Prime=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/3PrimeUTR
cupidPrediction=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred.3cols
cupidSeq=$cupidPrediction.sequence57bp
# cmd="$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/parseBS/getBSseq.py -s $seqof3Prime -p $cupidPrediction -o $cupidSeq &"
# $cmd

dnaseq=$CWD/refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv.uniq
# ~/bin/splitByN $dnaseq 1000
# for file in `ls ${dnaseq}_*`
# do
#   i=`echo $file|awk -F"/|_" '{print $NF}'`
#   cupidCoord=$file.cupidGenomicCoord
#   cmd="$PYTHON $srcDir/processData/parseBS/get3PrimeUTRCoords.py -c $cupidSeq  -d $file  -o $cupidCoord "
#   echo $cmd|qsub -l mem=16g,time=4:: -N bs_$i -cwd -e $CWD/log -o $CWD/log >> $CWD/qsub.log
#   tail -1 $CWD/qsub.log
# done

# awk 'NR==1||FNR!=1{print $0}' ${dnaseq}_*.cupidGenomicCoord > cupid.GenomicCoord

dnaseq=$CWD/tempfile/refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv.uniq
mycode="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step4-1_genMirBS.py"
utrseqfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/3PrimeUTR"
cupidfile="/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred"
output="/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/mircoRNA_BindSite_cupidPredict.hg19.test"
# for file in `ls ${dnaseq}_*`
# do
#   i=`echo $file|awk -F"/|_" '{print $NF}'`
#   cmd="$PYTHON $mycode -c $cupidfile  -u $utrseqfile -d $file  -o ${output}_$i"
#   # echo $cmd
#   echo $cmd|qsub -l mem=8g,time=4:: -N bs_$i -cwd -e $CWD/log -o $CWD/log >> $CWD/qsub.log
#   tail -1 $CWD/qsub.log
# done



