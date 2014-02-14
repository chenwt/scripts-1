#!/bin/bash
#$ -l mem=4g,time=20::
#$ -o ./log -e ./log
#$ -cwd
#Usage: runMuTect.sh <tumor.bam> <normal.bam> 
#Example
#Author:J.He
#Date: Jun 20,2013
##
echo "--------------This is the START-----------"
echo $(date)
bam1=$1'.rmdup.new.bam'
bam2=$2'.rmdup.new.bam'
cwd=$(dirname $0) 
bamDir='/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/callVars'
REF='/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa'
REFDBSNP='/ifs/scratch/c2b2/ac_lab/jh3283/ref/dbsnp_135_b37.vcf'
REFCOSMIC='/ifs/scratch/c2b2/ac_lab/jh3283/ref/b36_cosmic_v54_080711_sorted_2.vcf'
TUMOR=$bamDir'/'$bam1
NORMAL=$bamDir'/'$bam2
OUT=$cwd'/'$1_$2'_call_stats.txt'
#Intervals='/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/muTect/reTest_accepted.interval_list'
Intervals='/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/muTect/refSeq_codingExon_hg19.interval_list'
CovFile=$cwd'/'$1_$2'_coverage_wig.txt'
cmd="java -Xmx2g -jar /ifs/home/c2b2/ac_lab/jh3283/tools/muTect/muTect-1.1.4.jar
--analysis_type MuTect
--reference_sequence $REF
--dbsnp:VCF $REFDBSNP
--cosmic:VCF $REFCOSMIC
--intervals $Intervals
--input_file:normal $NORMAL
--input_file:tumor $TUMOR 
--out $OUT 
--vcf $OUT.vcf
--coverage_file $CovFile" 
echo $cmd
$cmd

##--------------------------------------
echo "start annotating....."
sh $ANNOVAR/do_annovar_all.sh $OUT.vcf $OUT.annotated.vcf


echo "----------This is the END--------"
echo $(date)
