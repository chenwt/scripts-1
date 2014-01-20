#!/bin/bash
#By: J.He
#TODO: 
#input:  coordinates to be search; bam file list; reference fa file path
#output: ATGC count of the input coordinates in the input bam files.

coorFile=$1
bamFile=$2
outputFile=$1"_"$2"_ACGT.txt"
echo -e "bamfile\tchr\tpos\tTotal\tTotal_Alt\tA\tC\tG\tT" > $outputFile
cnt_bam=0
for bam in `cat $bamFile`
do
  let "cnt_bam++"
  pid=`echo $bam|awk -F/ '{split($13,a,"-");print a[3]"-"a[4]}'` 
  echo "processing bamfile"$bam
  echo "/ifs/scratch/c2b2/ac_lab/rs3412/Scripts/novelSnvFilter_ACGT --var $coorFile --mapping $bam >temp"
  /ifs/scratch/c2b2/ac_lab/rs3412/Scripts/novelSnvFilter_ACGT --var $coorFile --mapping $bam >${pid}_ACGTrawCount.vcf
  #awk -v bam=$bam '{print bam"\t"$0}' temp >> $outputFile 
done

rm temp
echo "#------END------"
