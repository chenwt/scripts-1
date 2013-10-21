#!/bin/bash
#$ -cwd

workfolder=$1
bam=$2
ref=$3
targets=$4
outputfile=$5

maxdepth=10000
minMapQ=0
minBaseq=0


# the bed files have a 5 line header that must be removed prior to calling 
# samtools

targetstemp=${workfolder}/$(basename $targets).tmp.txt
tail -n +6 $targets > $targetstemp

cmd="samtools mpileup -DS -C 5  -d $maxdepth -q $minMapQ -Q $minBaseq  -f $ref -l $targetstemp $bam"

cmd=`$cmd > $outputfile`
$cmd
if [ $? -ne 0 ]; then
  echo "samtools did not exit properly"
  rm $targetstemp
  exit 1
fi
#rm $targetstemp
exit
