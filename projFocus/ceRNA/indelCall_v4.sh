#!/bin/bash
#By: J.He
#TODO: 1. ADD file existing checking for each step; 
#      2 ADD error handling for each step;  
#      3 ADD log file check for each step ;  
#      4 delete temp file after running. 
#Desp: created for projFocus ceRNA, to call indel using dindel using TCGA brca bam files
#     software version: dindel: dindel-1.01-linux-64bit, python version: Python 2.7.5
#     UNDER DEVELOPMENT
#input: <working directory full path> <BAM file full path> 
#ouput:
#COMMENT: require dindel installed
#updates: add filtering to dindel VCF, and output final filtered VCF compared to v2 which only output VCF

usage='usage: $0 <output dir> <full path to bam> 
      example: $0 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam'
#echo -e $usage

##--small func---
printSysTime() {
echo -e "step:\t"$1
echo $(date)|awk '{print "sysTime:\t"$2,$3,$4}'
}

wd=$1
bam=$2

###------example---start-----
#wd="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test"
#bam="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam"
###------example----end-----

##----initialization----
bamName=`echo $bam|awk 'BEGIN{FS="/"}{print $NF}'`
dindel="/ifs/home/c2b2/ac_lab/jh3283/tools/dindel/binaries/dindel-1.01-linux-64bit"
dindelDir="/ifs/home/c2b2/ac_lab/jh3283/tools/dindel/dindel-1.01-python"
ref="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa"
cd ${wd}

if [ ! -d ${wd}/log ]; then mkdir log ; fi
if [ ! -d ${wd}"/temp-"$bamName ]
then
  mkdir temp-$bamName 
else
  cd $wd/temp-$bamName
  rm -R *
  rm *
  cd $wd
fi

if [ ! -d ${wd}"/temp-"$bamName"/log" ]
then
  cd ${wd}"/temp-"$bamName
  mkdir log
  cd ${wd}
fi

wdFinal=$wd
outFileFinal=${wdFinal}"/"${bamName}"_dindel_ouput.variantCalls.VCF"
wd=$wd"/temp-"$bamName
log=$wd"/log"
outFile=${wd}"/"${bamName}"_dindel_ouput"

cd ${wd}
#--print parameter info
echo -e "working dir\t"$(pwd)
echo -e "BAM file:\t"$bam
echo -e "Output:\t"$outFile
echo -e "Final Output:\t"$outFileFinal
echo -e "log dir:\t"$log
echo -e "Reference:\t"$ref
echo -e "dindel version:\t"$dindel

##---step-1
printSysTime 1 
$dindel --analysis getCIGARindels --bamFile $bam \
--outputFile $outFile --ref $ref

#----step-2
printSysTime 2
${dindelDir}/makeWindows.py --inputVarFile $outFile.variants.txt \
--windowFilePrefix $outFile.realign_windows --numWindowsPerFile 1000


###----step-3##----need to do paralle calling for each window
printSysTime 3
cntJobSubmitted=0 
touch jobFinished.log
for windowFile in `ls $outFile.realign_window*txt|grep -v glf`
do
  tempScript=${wd}/windowfile$cntJobSubmitted.sh
  echo -e "realign using window file:\t $windowFile"

  echo -e '#!/bin/bash' >> $tempScript 
  cmd="$dindel --analysis indels --doDiploid --bamFile $bam --ref $ref \
  --varFile $windowFile \
  --libFile $outFile.libraries.txt \
  --outputFile $windowFile "
  echo $cmd >> $tempScript  
  echo 'echo "WindowFileFinished" >> jobFinished.log' >>$tempScript

  qsub -l mem=8g,time=125:: -N "win"$cntJobSubmitted -e $log/ -o $log/ -cwd $tempScript >> $log/qsub.log
  tail -1 $log/qsub.log 
  let cntJobSubmitted=$cntJobSubmitted+1
done

#step--3.5----job monitering
printSysTime 3.5
cntLoop=1
cntFinishedJob=$(awk 'END{print NR}' jobFinished.log)

while [[ $cntJobSubmitted -ne $cntFinishedJob ]]
do
  echo -e "Total Job: $cntJobSubmitted; Job finished: $cntFinishedJob "
  echo $(date)
  cntLoop=`echo "$cntLoop+1"|bc`
  sleep 10m
  echo "waiting $cntLoop loops..."
  cntFinishedJob=$(awk 'END{print NR}' jobFinished.log)
done

echo "Checking file number" 
numFileFinished=0
cntLoop=1
while [[ ! $cntJobSubmitted == $numFileFinished ]] 
do
   sleep 1m 
   numFileFinished=`ls $wd/*.glf.txt|wc -l `
   echo "total files: $cntJobSubmitted. finihsed files: $numfileFinished" 
   echo 'waiting for all glf.txt...'
   echo "$cntLoop+1"|bc
done

echo "All glf files generated, start calling...."
ls $wd/*glf.txt | tr "\t" "\n" > $wd/sample.dindel_stage2_outputfiles.txt


###----step-4-result
printSysTime 4
cntLoop=1
${dindelDir}/mergeOutputDiploid.py --inputFiles $wd/sample.dindel_stage2_outputfiles.txt \
--outputFile $outFile.variantCalls.VCF --ref $ref 
printSysTime 5

###----step-5-filtering
PYTHON=~/tools/python/Python_current/python
inputDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall
outputDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/indelFinal

###----
$PYTHON ~/scripts/projFocus/ceRNA/filterIndel.py -i $outFile.variantCalls.VCF -o $outFile.variantCalls.VCF.filtered.VCF

echo -e "copying final VCF..."
cp $outFile.variantCalls.VCF $wdFinal/
cp $outFile.variantCalls.VCF.filtered.VCF $wdFinal/

echo -e "removing temp file..."
cd $wdFinal
rm -R $wd
echo -e "output file:\t"$outFile.variantCalls.VCF.filtered.VCF
echo -e "#----END-----"


