##!/bin/bash
##$ -cwd
## ~/scripts/school/compGenomics/runMuTect/muTect_test.sh 

CWD='/ifs/scratch/c2b2/ys_lab/jh3283/school/compGenomic/LAML/'
REF='/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/wu_build36.fasta'
REFDBSNP='/ifs/scratch/c2b2/ys_lab/jh3283/ref/dbsnp_132_b36_resorted.vcf'
REFCOSMIC='/ifs/scratch/c2b2/ys_lab/jh3283/ref/b36_cosmic_v54_080711_sorted_2.vcf'
TUMOR=$CWD'WGS/1bb5b7ef-16da-43ab-a8f4-2d98d1770e21/TCGA-AB-2978-03A-01D-0739-09_whole.bam'
NORMAL=$CWD'WGS/05f440e8-6e9d-43b5-9340-5271aed310dd/TCGA-AB-2972-11A-01D-0739-09_whole.bam'

cmd="java -Xmx2g -jar /ifs/home/c2b2/ys_lab/jh3283/tools/muTect/muTect-1.1.4.jar
--analysis_type MuTect
--reference_sequence $REF
--dbsnp $REFDBSNP
--cosmic $REFCOSMIC
--intervals 1:75100-7577200
--input_file:normal $NORMAL
--input_file:tumor $TUMOR 
--out test_call_stats.txt 
--coverage_file test_coverage_wig.txt"
echo $cmd
$cmd

##------------------------------test
java -Xmx2g -jar /ifs/home/c2b2/ys_lab/jh3283/tools/muTect/muTect-1.1.4.jar
--analysis_type MuTect
--reference_sequence REF
--dbsnp dbsnptest.vcf
--cosmic COSMICVCF
--intervals 1:75100-7577200
--input_file:normal normal.bam
--input_file:tumor tumor.bam 
--out test_call_stats.txt 
--coverage_file test_coverage_wig.txt