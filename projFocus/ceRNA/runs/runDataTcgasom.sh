#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 

###test
# head -100 genome.wustl.edu__Illumina_All.maf.matrix > toyfile
# cut -f1 toyfile|sort|uniq >toygene
# regfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_coords_Tss.tsv.signle.tsv
# toyreg=toyfile.region.tss
# echo -n "" > $toyreg
# while read gene
# do
#   awk -v g=$gene '$1==g{print $0}' $regfile >> $toyreg
# done < toygene

# sort $toyreg |uniq > $toyreg.region.tss.single
# input=toyfile
# tarregf=toyfile.region.tss.region.tss.single
# promotersize=2000
# output=toyout
# ~/tools/python/Python_current/python $crnsc/processData/selectTargetRegionRow.py -i $input -t $tarregf -c $promotersize -o $output
# rm toy*


CWD='/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som'
CDT=`date|awk '{print $2"-"$3"-"$6}'`
input=$CWD/genome.wustl.edu__Illumina_All.maf.matrix
# tssfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_coords_Tss.tsv.signle.tsv
tssfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar20_Tss.tsv.signle.tsv
promotersize=2000
output=$input.promoter2k.$CDT.matrix
# cmd="~/tools/python/Python_current/python $crnsc/processData/selectTargetRegionRow.py -i $input -t $tssfile -c $promotersize -o $output &"
# echo $cmd
# ~/tools/python/Python_current/python $crnsc/processData/selectTargetRegionRow.py -i $input -t $tssfile -c $promotersize -o $output &
# promotersize=1000
# output=$input.promoter1k.$CDT.matrix
# ~/tools/python/Python_current/python $crnsc/processData/selectTargetRegionRow.py -i $input -t $tssfile -c $promotersize -o $output &

# utr3pfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_coords_3pUTR.tsv
utr3pfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar20_3pUTR.tsv
size=0
output=$input.utr3p.$CDT.matrix
# ~/tools/python/Python_current/python $crnsc/processData/selectTargetRegionRow.py -i $input -t $utr3pfile -c $size -o $output &


# utr5pfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_coords_5pUTR.tsv
utr5pfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar20_5pUTR.tsv
size=0
output=$input.utr5p.$CDT.matrix
# ~/tools/python/Python_current/python $crnsc/processData/selectTargetRegionRow.py -i $input -t $utr5pfile -c $size -o $output &


##---getRecurrence 
# ~/scripts/projFocus/ceRNA/plot/getMutFreqInSamples.sh genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix &
## sort -k 5n genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix.sorted.mutFreq|tac |less
#sortxh /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix > /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-20-2014.matrix.sorted
#/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/plot/getMutFreqForGene.sh brca_som_utr3p.Mar-20-2014.matrix.mutFreq.sorted
