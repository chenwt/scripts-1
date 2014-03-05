#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: run breakdancer v1.1r81

bdDir=/ifs/home/c2b2/ac_lab/jh3283/tools/breakdancer/breakdancer-1.1.2
cwd=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/chrRecomb/
bamDir=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/callVars/
tumorBam="PARGVC-03A.rmdup.new.bam" 
relapseBam="PARGVC-04A.rmdup.new.bam"
normalBam="PARGVC-14A.rmdup.new.bam"

pid="PARGVC"
if [[ ! -d $cwd/temp-$pid ]]; then mkdir $cwd/temp-$pid ; fi
if [[ ! -d $cwd/temp-$pid/log ]]; then mkdir $cwd/temp-$pid/log ; fi

tempDir=$cwd/temp-$pid
cd $tempDir
$bdDir/perl/bam2cfg.pl -g -h $bamDir/$tumorBam $bamDir/$normalBam > $tempDir/${pid}_tn.cfg
perl /ifs/home/c2b2/ac_lab/jh3283/tools/breakdancer/breakdancer/branches/release-r76/BreakDancerMax.pl -g ${pid}_tn.bed -t -q 10 -d ${pid}_tn.ctx ${pid}_tn.cfg > ${pid}_tn.ctx

