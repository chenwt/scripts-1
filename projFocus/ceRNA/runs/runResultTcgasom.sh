#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

PYTHON=~/tools/python/Python_current/python
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/tcgal2som
# cut -f2 genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix.groupbyGene.freq|sort|uniq -c

# cp genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix.sorted brca_tcga_som.promoter2k.sorted.matrix 
# cp genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-20-2014.matrix brca_tcga_som.utr3p.matrix

utr3mat=$CWD/brca_tcga_som.utr3p.matrix
# cut -f 1,2,6- $utr3mat |sortxh |uniq > $utr3mat.sorted.uniq
utr3mat=$CWD/brca_tcga_som.utr3p.matrix.sorted.uniq

promotermat=$CWD/brca_tcga_som.promoter2k.sorted.matrix
# cut -f 1,2,6- $promotermat |sortxh |uniq > $promotermat.sorted.uniq
promotermat=$CWD/brca_tcga_som.promoter2k.sorted.matrix.sorted.uniq

output1=$utr3mat.groupbyGene.freq
output2=$promotermat.groupbyGene.freq
# # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/tcgaSom/groupMut.py -i $utr3mat -o $output1  
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/tcgaSom/groupMut.py -i $promotermat -o  $output2 

# sort -k 2nr $output2 |cut -f2 |sort |uniq -c |sort -k 2n > $output2.table
sort -k 2nr $output1 |cut -f2 |sort |uniq -c |sort -k 2n > $output1.table

## get union of promoter and 3'UTR
# awk 'NR==1||FNR!=1{print $0}' brca_tcga_som*matrix  > brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix 
# ~/bin/sortxh brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix|uniq > brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted
#cut -f 1,2,6- brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted |sortxh |uniq > brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.uniq
allmat=$CWD/brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.uniq
output3=$allmat.groupbyGene.freq
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/tcgaSom/groupMut.py -i $allmat -o  $output3 
# sort -k 2nr brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.groupbyGene.freq |cut -f2 |sort |uniq -c |sort -k 2n > brca_tcga_som_promotor_utr3p_Apr-9-2014.matrix.sorted.groupbyGene.freq.table

# sort -k 2nr $output3 |cut -f2 |sort |uniq -c |sort -k 2n > $output3.table
