#!/bin/bash
#By: J.He
#TODO: 
##--------functions-------------
genInputCall(){
  dataDir=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/
  input=$1
  output=$2
  while read line
  do
    barcode=`echo $line|awk '{print $1}'`
    readlink -f $dataDir/$barcode.bam >> $output
  done < $input
}


function doCallVars(){
  cwd=`pwd`
  if [ ! -d $cwd/log ]; then mkdir $cwd/log ; fi
  logDir=$cwd/log
  cnt=0
  while read bam 
  do
    ((cnt+=1))
    pid=`echo $bam|awk -F"/" '{split($NF,a,"-");print a[3]}'`
    echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-1_splitBAM.sh $bam" |qsub -l mem=2g,time=140:: -N fireAll_$pid -cwd -e $logDir -o $logDir >> qsub.log
    tail -1 qsub.log
  done < $1
  echo "$cnt sample start calling"
}

#----------------executions-------------
# cp /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part1 .
# cp /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part11 . 
# genInputCall input_gtBatch_v2.txt_part11 input_wgsCall_part11.txt
# head -1 input_wgsCall_part11.txt> input_wgsCall_partTest.txt
# sed -i 1d input_wgsCall_part11.txt 
# head -1 input_wgsCall_part11.txt> input_wgsCall_part0.txt
# doCallVars input_wgsCall_partTest.txt

# doCallVars input_wgsCall_part0.txt

echo "##---END---"
