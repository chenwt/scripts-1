#!/bin/bash
#By: J.He
#$ -cwd 
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
  echo -e "$input" > $input.done
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
    # echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-1_splitBAM.sh $bam" |qsub -l mem=2g,time=140:: -N fireAll_$pid -cwd -e $logDir -o $logDir >> qsub.log
    # tail -1 qsub.log
    /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/doSubBAMcall.sh $bam
    # sleep 60m
    sleep 120m
  done < $1
  echo "$cnt sample start calling"
}

function mergeAll(){
  ###-----------do merging ------
  tempDir_array=( `ls -r |grep "temp" ` )
  # echo ${tempDir_array[@]}
  for tempDir in "${tempDir_array[@]}"
  do
      # echo $tempDir
      cntVCF=`ls $tempDir/*var.vcf.gatk.vcf 2>/dev/null |wc -l`
      cntBAM=`ls $tempDir/split*.bam  2>/dev/null |wc -l`
      # cntTotal=`wc -l $REGIONFILE 2>/dev/null|wc -l `
      cntTotal=103
      echo "$tempDir vcf $cntVCF bam $cntBAM total $cntTotal" 
      if [ $cntVCF != $cntTotal ] || [ $cntBAM -gt 0 ]
      then
          # echo -e "skip merging $tempDir "
	  continue 
      fi
      echo "merging...." 
      pid=`echo $tempDir|awk -F"-" '{print $4}' `
      dir=`readlink -m $tempDir`
      # echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/ " |qsub -l mem=8g,time=4:: -N merge_$pid -e ./log -o ./log -cwd 
      
      echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeAnnoFilter.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/ " |qsub -l mem=12g,time=8:: -N merge_$pid -e ./log -o ./log -cwd 
      # /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeAnnoFilter.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/ & 
  done 
}

delBAMCalled(){
  ##------deleting bam files-------
  resDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/ 
  dataDir=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS
  /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/delBAM.sh $resDir $dataDir
}

#----------------executions-------------
###-------DATA
 # delBAMCalled 
 mergeAll

# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part10 input_wgsCall_part10.txt
# doCallVars input_wgsCall_part10.txt &

# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part9 input_wgsCall_part9.txt
# doCallVars input_wgsCall_part9.txt &

 # genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part8_section2Feb25 input_wgsCall_part8_section2.txt
 # doCallVars input_wgsCall_part8_section2.txt &
# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part8_section1Feb25 input_wgsCall_part8_section1.txt
 # doCallVars input_wgsCall_part8_section1.txt & 
 # doCallVars input_wgsCall_part6_rescue_Feb242014.txt 
# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part7  input_wgsCall_part7.txt
# doCallVars input_wgsCall_part7.txt &
# doCallVars input_wgsCall_rescue_Feb232014.txt &
# cp /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part1 .
# cp /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part11 . 
# genInputCall input_gtBatch_v2.txt_part11 input_wgsCall_part11.txt
# head -1 input_wgsCall_part11.txt> input_wgsCall_partTest.txt
# sed -i 1d input_wgsCall_part11.txt 
# head -1 input_wgsCall_part11.txt> input_wgsCall_part0.txt
# sed -i 1d input_wgsCall_part11.txt 
# head -1 input_wgsCall_part11.txt >  input_wgsCall_part00.txt
# sed -i 1d input_wgsCall_part11.txt 

# doCallVars input_wgsCall_partTest.txt
# doCallVars input_wgsCall_part0.txt
# doCallVars input_wgsCall_part00.txt
# doCallVars input_wgsCall_part11.txt
# doCallVars input_wgsCall_part2.txt
# doCallVars input_wgsCall_part_normal_3.txt 

#genInputCall input_gtBatch_v2.txt_part11 input_wgsCall_part1.txt
# doCallVars input_wgsCall_part11.txt

# genInputCall input_gtBatch_v2.txt_part1 input_wgsCall_part1.txt
# doCallVars input_wgsCall_part1.txt

# doCallVars input_wgsCall_part10.txt

# doCallVars input_wgsCall_part1.txt
# doCallVars input_wgsCall_part1_A04Q.txt
# doCallVars input_wgsCall_part_normal_3.1.txt 

# cp /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part2 . 
# genInputCall input_gtBatch_v2.txt_part2 input_wgsCall_part2.txt

# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part3 input_wgsCall_part3.txt
# doCallVars input_wgsCall_part3.txt &

# mv input_wgsCall_part3.txt done/
# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part5 input_wgsCall_part5.txt
# doCallVars input_wgsCall_part5.txt &

####-----TCGA----
# tcgaCallDir=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/CALL
# cd $tcgaCallDir
# ln -s /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/test_doSteps.sh RUN.sh
# cp $fdata/wgs/input_wgsCall_part_normal_5.txt .
# doCallVars input_wgsCall_part_normal_5.txt 
# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part2 input_wgsCall_part2.txt
# doCallVars input_wgsCall_part21.txt 
# doCallVars input_wgsCall_part2.txt 
# doCallVars input_wgsCall_part2_rsq.txt 
# doCallVars input_wgsCall_part_normal_5.1.txt 
 # sleep 100m
 # doCallVars input_wgsCall_part_normal_5.txt &
# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part4 input_wgsCall_part4.txt
# doCallVars input_wgsCall_part4.txt  &
# doCallVars input_wgsCall_part4_A0D0.txt & 

# genInputCall /ifs/scratch/c2b2/TCGA/data/BRCA/WGS/input_gtBatch_v2.txt_part2.redownload input_wgsCall_part2.redownload.txt
# doCallVars input_wgsCall_part2.redownload.txt  &

 # mergeAll

echo "##---END---"
