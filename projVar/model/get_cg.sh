#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

mfile=$1
bin=$2 ## bin in 10, 100, 1000, 10000, 100000
ofile=$mfile.$bin.cg

$BTDIR/bedtools intersect -wb -a -b  

