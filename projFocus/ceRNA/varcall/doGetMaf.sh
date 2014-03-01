#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test
cwd="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test"
# head -500 ../../rawVars/TCGA-BH-A0E0-11A-13D-A128-09.bam.var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf > test.vcf
PYTHON=~/tools/python/Python_current/python

vcfdir='/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test'
# genelist='/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/cancerGenes_UCceRNETPlusNature12929.list.geneSingleStartEnd'
# genelist='test.gene.list'
# keyType='maf' 
# locType='tss'
# example='python " + sys.argv[0] + " -d $vcfdir -o $output -l $locType -k $keyType -g $genelist'
# output="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test/test.${locType}.${keyType}.out"
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-6_getMAF.py  -d $cwd -o $output -l $locType -k $keyType -g $genelist  

genelist='test.gene.list'
locType='gene'
keyType='snp'
output="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/mafs/test/test.${locType}.${keyType}.out"
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-6_getMAF.py  -d $cwd -o $output -l $locType -k $keyType -g $genelist  
