#!/bin/bash
#By: J.He
#TODO: 
##input: downloaded file name(from input_gtDownload) ; tcga tsv information file
##output: problematic files which haven't finished 
##example; $0 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/wgs/input_gtBatch_v2.txt_part2 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/wgs/brca_wgs_bam_summary_02042014.tsv   
dataDir=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS

echo -n "" > $1.redownload
while read line
do 
  barcode=`echo $line|awk '{print $1}'`
  tempCntFile=`ls $dataDir/$barcode.bam 2>/dev/null|wc -l`
  if [[ $tempCntFile == 0 ]]
  then
    echo "$barcode.bam not exisit"
    echo $line >> $1.redownload
    continue
  else
    size=`grep $barcode $2|awk -F"\t" '{print $15}'`
    dSize=`ls -l $dataDir/$barcode.bam|awk '{print $5}'`
    if [[ $dSize != $size ]]
    then
      echo "download $dSize actually size $size $barcode"
      echo $line >> $1.redownload
      continue
    fi
    echo $barcode.bam " download success"  

  fi
done < $1
