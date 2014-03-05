#!/bin/bash
#! -cwd
#By: J.He
#Desp.: misc operations on  

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth
awk 'FNR!=1||NR==1' $CWD/brca_methTumor_level3_02072014.mat_diffMeth.mat \
$CWD/brca_methTumor_nature12912NovelGene_02172014.mat_diffMeth.mat \
|~/bin/sortxh |uniq> $CWD/brca_methTumor_combinedCG_03032014.mat_diffMeth.mat

