#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
CWD=/ifs/home/c2b2/ac_lab/jh3283/DATA/projMisc/varReg/data/062014/tcgaSom
##step 1 positive control
# python /ifs/home/c2b2/ac_lab/jh3283/scripts/projVar/processData/step1_extractMut.py -b ../wgsVar/brca_wgs_paired_tumor.barcode.list -m genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf  -o brca_tcga_70tumor_wustlexmut.postive.mutation

posMut=$CWD/brca_tcga_70tumor_wustlexmut.postive.mutation
awk 'BEGIN{OFS="\t"}{print "chr"$1,$2,$3,$4,$5}' $posMut > $CWD/positive.mutation


##step 2 negative control
posMut=$CWD/positive.mutation
negMut=$CWD/negative.mutation
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projVar/model/step1_genRandMut.py -s $posMut -o $negMut


