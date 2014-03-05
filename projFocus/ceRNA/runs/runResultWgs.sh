#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: scratches for running jobs

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered
in=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/wgs/input_gtBatch_v2_WU.txt

cd $CWD
out=WU.vcflist
echo -n "" > WU.vcflist
while read line 
do
  temp=`echo $line|awk '{print $1}'`
  cnt=`ls $temp*vcf 2>/dev/null|wc -l` 
  if [ $cnt -gt 0  ]; then  
    ls -1 $temp*vcf >> $out
    mv $temp*vcf $CWD/wu/
  fi
done < $in

