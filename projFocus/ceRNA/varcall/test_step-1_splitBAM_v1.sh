#!/bin/bash
#By: J.He
#$ -cwd 
#TODO: 


##----split bam file
GATKJAR=/ifs/home/c2b2/ac_lab/jh3283/tools/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
# REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/humanGenome_hg19_regions104.bed
REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/test/hg19_region_test.bed
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
BAM=/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/TCGA-A8-A08B-01A-11D-A19H-09.bam
bamName=`echo $BAM|awk 'BEGIN{FS="/"}{print $NF}'`
if [ ! -d temp-$bamName ]; then mkdir temp-$bamName ; fi
rootd=`pwd`
cwd=$rootd/temp-$bamName
if [ ! -d $cwd/log ]; then mkdir $cwd/log ; fi
logDir=$cwd/log
###--need 3-4 hours for one 
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/splitByRegion.sh $REGIONFILE $BAM $cwd

####----------
# #invoke call variants job
cd $cwd

# sleep 30m 
# sleep 30s 
cntCallJob=0
numSplitJobDone=0
numSubBAM=`awk 'END{print NR}' $REGIONFILE`
echo -e "debug numSubBAM $numSubBAM"
if [ ! -d $cwd/callJobSubmitted ]; then mkdir $cwd/callJobSubmitted ; fi 
while [ $cntCallJob -le $numSubBAM ]; do
  numSplitJobDone=`ls $cwd/split_*bam.list|awk 'BEGIN{FS="\t"}END{print NR}'`
  echo "debug $numSplitJobDone"
  if [[ $numSplitJobDone -gt 0 && $cntCallJob -le $numSubBAM ]] 
  then
     for list in $( ls $cwd/split_*bam.list); do
        subBam=`echo $list|awk '{gsub(".list","",$1);print $0}'`
	echo "debug $subBam"
	subBam=`readlink -f $subBam`
	echo "debug bam: $subBam"
	cmd="$srcDir/test_step-2_callVar.sh $cwd $subBam "
	echo $cmd
        mv $list $cwd/callJobSubmitted/ 
	let cntCallJob=$cntCallJob+1
	echo "debug calljob: $cntCallJob"
    done
    echo "debug $cntCallJob"
    # exit
  else
    echo -e "no split*bam.list found"
  fi
  # sleep 30m
  sleep 30s
done

if [  ! -f $cwd/samtools_callSomaticVar.sh ] ;then  
  cp $srcDir/samtools_callSomaticVar.sh $cwd 
    chmod +x $cwd/samtools_callSomaticVar.sh
fi

  for cnt in `seq 1 $numSubBAM` 
  do
     cmd="qsub -N myjobCall$cnt -e $logDir -o $logDir -b y -l mem=8G,time=8:: -cwd \"$cwd/samtools_callSomaticVar.sh $cwd/split_$cnt.$bamName > jobFinished.txt\" "
##     $cmd > $cwd/qsub.log
##     tail -1 $cwd/qsub.log
##  done
