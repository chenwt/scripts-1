#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

##step 1 positive control
# python /ifs/home/c2b2/ac_lab/jh3283/scripts/projVar/processData/step1_extractMut.py -b ../wgsVar/brca_wgs_paired_tumor.barcode.list -m genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf  -o brca_tcga_70tumor_wustlexmut.postive.mutation


##step 2 negative control
posMut=$CWD/brca_tcga_70tumor_wustlexmut.postive.mutation
negMut=$CWD/negative.mutation
python /ifs/home/c2b2/ac_lab/jh3283/scripts/projVar/model/step1_genRandMut.py -s $posMut -o $negMut


