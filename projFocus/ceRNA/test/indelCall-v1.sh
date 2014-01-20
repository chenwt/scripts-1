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

usage='usage: $0 <output dir> <full path to bam> 
      example: $0 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam'
#echo -e $usage

##--small func---
printSysTime() {
echo "step:\t"$1
echo $(date)|awk '{print "sysTime:\t"$2,$3,$4}'
}

##------example---start-----
wd="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test"
bam="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam"
##------example----end-----

##----initialization----
bamName=`echo $bam|awk 'BEGIN{FS="/"}{print $NF}'`
dindel="/ifs/home/c2b2/ac_lab/jh3283/tools/dindel/binaries/dindel-1.01-linux-64bit"
dindelDir="/ifs/home/c2b2/ac_lab/jh3283/tools/dindel/dindel-1.01-python"
ref="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa"
cd ${wd}
if [ ! -d ${wd}/log ]; then mkdir log ; fi
if [ ! -d ${wd}/temp ]; then mkdir temp ; fi

wdFinal=$wd
outFileFinal=${wd}"/"${bamName}"_dindel_ouput.variantCalls.VCF"
log=$wc"/log"
wd=$wd"/temp"
outFile=${wd}"/"${bamName}"_dindel_ouput"

#--print parameter info
echo -e "working dir\t"$wd
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

#####----step-2.5-select_conditional_candidates if necessary
####---this sub step was not included in routine calling
##printSysTime 2.5
##${dindelDir}/selectCandidates.py -i sample.dindel_output.variants.txt \
##-o sample.dindel_output.variants.mincount2.txt


###----step-3##----need to do paralle calling for each window
printSysTime 3
cntJobSubmitted=0 
for windowFile in `ls $outFile.realign_window*txt|grep -v glf `
do
  echo -e "realign using window file:\t $windowFile"
  cmd="$dindel --analysis indels --doDiploid --bamFile $bam --ref $ref \
  --varFile $windowFile \
  --libFile $outFile.libraries.txt \
  --outputFile $windowFile "
  #jobName=`echo $outFile|awk 'BEGIN{FS="\\"}{print $NF}'`
  echo $cmd |qsub -l mem=8g,time=24:: -N "windowFile"$cntJobSubmitted -e $log/ -o $log/ -cwd  
  cntJobSubmitted=`echo "$cntJobSubmitted+1"|bc`
done

##windowFile=$wd/test.bam_dindel_ouptu.realign_windows.10.txt
windowFile=`ls $outFile.realign_window*txt|grep -v glf|tail -1`
#step--3.5----job monitering
printSysTime 3.5
cntLoop=1
while [ ! -e $windowFile".glf.txt" ] 
do
  echo "waiting for job...."
  echo $(date)
  echo "$cntLoop+1"|bc
  sleep 5m
  echo "waiting 5 * $cntLoop mins"
done

echo "here" 
numFileFinished=`ls $wd/*.glf.txt|wc -l `
cntLoop=1
while [[ ! $cntJobSubmitted == $numFileFinished ]] 
do
  echo "total files: $cntJobSubmitted. finihsed files: $numfileFinished" 
   echo 'waiting for all glf.txt...'
   echo "$cntLoop+1"|bc
   sleep 1m ##
done

echo "All glf files generated, start calling...."
ls $wd/*glf.txt | tr "\t" "\n" > $wd/sample.dindel_stage2_outputfiles.txt


###----step-4-result
printSysTime 4
cntLoop=1
${dindelDir}/mergeOutputDiploid.py --inputFiles $wd/sample.dindel_stage2_outputfiles.txt \
--outputFile $outFile.variantCalls.VCF --ref $ref 
printSysTime END
echo -e "copying final VCF..."
cp $outFile.variantCalls.VCF $wdFinal/
echo -e "removing temp file..."
rm -R $wd
echo -e "output file:\t"$outFile.variantCalls.VCF
echo -e "#----END-----"
