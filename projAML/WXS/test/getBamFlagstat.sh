#!/bin/bash
#By: J.He
#TODO: 
#input:<file:list of full path of bam files to get flagstat>
#output:<file:flag stat for each bamfile>
samtool=`which samtools`
outfile=$1".flagstat.txt"
if [ -f $outfile ]; then rm $outfile ; fi
cnt_file=0
num_file=`wc -l <$1`
for bam in `cat $1`
do
  let cnt_file++
  echo "processing file "$cnt_file" : "$bam"..."
  if [ $cnt_file==$num_file ]; then
    echo "last file"
    $samtool flagstat $bam > $cnt_file.temp 
  else 
    $samtool flagstat $bam > $cnt_file.temp &
  fi
done

echo -e "in_total\tduplicates\tmapped\tpaired\tread1\tread2\tproperly_paired\twith_itself_and_mate_mapped\tsingletons\twith_mate_mapped_diff_chr\twith_mate_mapped_diff_chr_mapQ5" > $outfile
for cnt in `seq 1 $cnt_file`
do
  cat $cnt.temp | awk 'BEGIN{ORS="\t"}{print $1" + "$3}'|sed 's/\t$/\n/'  >> $outfile
done


echo "#-------------END-------"
