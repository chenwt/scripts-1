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
~/bin/splitByN $dnaseq 3000
# cupidCoord=$cupidPrediction.genomicCoord_Apr6
# $PYTHON $srcDir/processData/parseBS/get3PrimeUTRCoords.py -c $cupidSeq  -d $dnaseq  -o $cupidCoord 
