#!/bin/bash
#$ -cwd
#By: J.He
#Desp.: 


# awk '{print $1, $2, $4}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/CupidPred | uniq > CupidPred.3cols

/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/parseBS/fasta2tsv.py refseq_hg19_refflat_cupid3pGene.fasta refseq_hg19_refflat_cupid3pGene.fasta_Apr4.tsv 
