#!/bin/bash
#By: J.He
#TODO: 
#Desp: created for projFocus ceRNA, for call indel using dindel using TCGA brca bam files
#     UNDER DEVELOPMENT
#input:
#ouput:
#COMMENT: require dindel installed

##----initialization----
wd="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test"
dindel="/ifs/home/c2b2/ac_lab/jh3283/tools/dindel/binaries/dindel-1.01-linux-64bit"
dindelDir="/ifs/home/c2b2/ac_lab/jh3283/tools/dindel/dindel-1.01-python"
bam="/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam"
ref="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa"
outFile=$bam"_dindel_ouptu"
####---step-1
##$dindel --analysis getCIGARindels --bamFile $bam \
##--outputFile $outFile --ref $ref
#

###----step-2
##${dindelDir}/makeWindows.py --inputVarFile $outFile.variants.txt \
##--windowFilePrefix $outFile.realign_windows --numWindowsPerFile 1000

###----step-2.5-select_conditional_candidates if necessary
#${dindelDir}/convertVCFToDindel.py  
#${dindelDir}/selectCandidates.py -i sample.dindel_output.variants.txt \
#-o sample.dindel_output.variants.mincount2.txt###


##----step-3##----need to do paralle calling for each window
cntJobSubmitted=0 
for windowFile in `ls $outFile.realign_window*txt|grep -v glf |head -5`
do
  echo -e "realign using window file:\t $windowFile"
  cmd="$dindel --analysis indels --doDiploid --bamFile $bam --ref $ref \
  --varFile $windowFile \
  --libFile $outFile.libraries.txt \
  --outputFile $windowFile "
  #jobName=`echo $outFile|awk 'BEGIN{FS="\\"}{print $NF}'`
  echo $cmd |qsub -l mem=8g,time=24:: -N "windowFile"$cntJobSubmitted -cwd  
  cntJobSubmitted=`echo "$cntJobSubmitted+1"|bc`
done

#step--3.5----job monitering

while [ ! -e $outFile ] 
do
  sleep 30m

done

numFileFinished=`ls $wd/*.glf.txt|wc -l `
while [[ ! $cntJobSubmitted == $numFileFinished ]] 
do
   sleep 1m ##
   echo 'waiting for glf.txt'
done

echo "All glf files generated, start calling...."
ls $wd/*glf.txt | tr "\t" "\n" > $wd/sample.dindel_stage2_outputfiles.txt

###----step-4-result
${dindelDir}/mergeOutputDiploid.py --inputFiles $wd/sample.dindel_stage2_outputfiles.txt \
--outputFile $outFile.variantCalls.VCF --ref $ref 

