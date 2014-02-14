#!/bin/bash
#By: J.He
#$ -cwd 
#TODO: 
##input: <full path of bam > <bed file include the region to split the bam>
##output:

BAM=$1
##----split bam file
GLOBAL=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh
. $GLOBAL
GATKJAR=/ifs/home/c2b2/ac_lab/jh3283/tools/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/humanGenome_hg19_regions104.bed
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
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
if [ ! -d $cwd/callJobSubmitted ]; then mkdir $cwd/callJobSubmitted ; fi 
cnt=`awk 'END{print NR}' $REGIONFILE`
((cnt -= 1))
echo -e "total region file:\t"$cnt
cntNF=`ls -1 split_*bam.list 2>/dev/null | wc -l`
cntjob=0
echo "new file.." $cntNF
while [ $cntNF -lt $cnt ] 
      [ $cntjob -lt $cnt ]
do
    if [ $cntNF -gt 1 ] ; then
      for list in `ls split_*bam.list` 
      do
         subBam=`readlink -f $list|awk '{gsub(".list","",$1);print $1}'`
	 cmd="$srcDir/step-2_callVar.sh $cwd $subBam"
	 echo $cmd | qsub -l mem=8g,time=24:: -N step2_$cntjob -e $logDir -o $logDir -cwd  >>$cwd/qsubStep2.log
         mv $list $cwd/callJobSubmitted/ 
	 ((cntjob += 1))
         echo "$list calling started" 
	 echo "$cntjob out of $cnt split bam done.."
      done
    elif [ $cntNF -eq 1 ]; then 
      list=`ls split_*bam.list`
      subBam=`readlink -f $list|awk '{gsub(".list","",$1);print $1}'`
      cmd="$srcDir/step-2_callVar.sh $cwd $subBam"
      echo $cmd | qsub -l mem=8g,time=24:: -N step2_$cntjob -e $logDir -o $logDir -cwd >>$cwd/qsubStep2.log
      mv $list $cwd/callJobSubmitted/ 
      ((cntjob += 1))
      echo "$list calling started" 
      echo "$cntjob out of $cnt split bam done.."
    else
      echo "$cntjob out of $cnt split bam done....$(date)"
    fi
    if [ $cnt -eq $cntjob ]; then break ; fi
    sleep 10m
    cntNF=`ls -1 split_*bam.list 2>/dev/null | wc -l`
done

echo "All calling job submitted!"
echo "#---DONE---split-and-start-calling-----"
