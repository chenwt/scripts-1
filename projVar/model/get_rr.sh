#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

mfile=/ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/062014/tcgaSom/positive.mutation
# mfile=$1

rrfile=/ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/genomicFeature/rr/hapMapRelease24CombinedRecombMap.bw.wig
rrXfile=/ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/genomicFeature/rr/Female.bw.wig

BTDIR=/ifs/home/c2b2/ac_lab/jh3283/tools/bedtools-2.17.0/bin/

ofile=$mfile.rr.bed

$BTDIR/bedtools intersect -wb -a $mfile -b $rrfile > temp 
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"rr",$9}' temp > $ofile 

rm temp

$BTDIR/bedtools intersect -wb -a $mfile -b $rrXfile > temp 
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,"rr", $9}' temp >> $ofile 

# rm temp
echo "[END]"

