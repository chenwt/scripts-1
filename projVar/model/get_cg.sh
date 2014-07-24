#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

BTDIR=/ifs/home/c2b2/ac_lab/jh3283/tools/bedtools-2.17.0/bin/
mfile=$1

for bin in 10 100 1000 10000 100000
do
  cgfile=/ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/genomicFeature/gc/gcPercent_${bin}.bed
  ofile=$mfile.cg.$bin.bed
  # echo $bin
  # echo $cgfile
  $BTDIR/bedtools intersect -wb -a $mfile -b $cgfile > temp 
  awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$9, $10}' temp > $ofile 
done 
rm temp
echo "[END]"

