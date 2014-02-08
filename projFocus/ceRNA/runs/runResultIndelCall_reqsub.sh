#!/bin/bash
#By: J.He
#TODO: 

###----test-----

#echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v2.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam"|qsub -l mem=8g,time=40:: -N indelTestrun -cwd -o ./log -e ./log >>qsub.log 
#tail -1 qsub.log

#data folder:
dataDir=/ifs/scratch/c2b2/TCGA/data/BRCA/WXS
outDir=$(pwd)
##-----test for one sample
##TCGA-A1-A0SD-01A-11D-A10Y-09.bam
#echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v2.sh $outDir /ifs/scratch/c2b2/TCGA/data/BRCA/WXS/TCGA-A1-A0SD-01A-11D-A10Y-09.bam "|qsub -l mem=8g,time=40:: -N indelTestrun -cwd -o ./log -e ./log >>qsub.log
#tail -1 qsub.log

###---run for listed files
## at  most 20 per run
##--making bamlist files
#cntTotalBam=`ls ${dataDir}/*bam|awk 'END{print NR}'`
#cnt=`echo "$cntTotalBam/20+1"|bc`
#cntPart=1
#echo -e "$cntTotalBam\t$cnt\t$cntPart"
#ls ${dataDir}/*bam > inputIndelCallBamList 
#for i in `seq 1 $cnt`
#do
# cntx=`echo $i*20|bc`
# cntm=`echo "$i*20-20"|bc`
# echo -e "line start:$cntm\tline end:$cntx"
# awk -v cntx=$cntx -v cntm=$cntm 'NR<=cntx&&NR>cntm{print $0}' inputIndelCallBamList > input.bam.list.${i}
#done

#submit run
qsubRun (){
  cnt=0
  while read line
  do
    let cnt=$cnt+1
    jobid=`echo $line|awk 'BEGIN{FS="/"}{print substr($NF,9,4)}'`"_$cnt"
    echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3.sh $outDir $line"|qsub -l mem=8g,time=140:: -N $jobid -cwd -o ./log -e ./log >>qsub.log
    tail -1 qsub.log
  done <$1
}


###---first run, time only 40 hrs
#qsubRun input.bam.list.7 
#qsubRun input.bam.list.0 
#qsubRun input.bam.list.1 

###debugging for first run
###waiting for all windown file job finished 
#for line in `cat input.bam.list.1 input.bam.list.0`
check() {
  while read line
  do
    bam=`echo $line|awk 'BEGIN{FS="/"}{print $NF}'` 
    if [ -d temp-$bam ] ;then
      cd temp-$bam
      jobFinished=`awk 'END{print NR}' jobFinished.log `
      jobSubmitted=`ls win*sh|awk 'END{print NR}'`
      if [ $jobFinished -eq $jobSubmitted ] ; then
        echo "$bam finishined..."
        echo $bam >> ../bam.finished.txt
      else 
        jobRun=`echo $jobSubmitted-$jobFinished|bc`
        echo -e "$bam \t job running $jobRun"
      fi
      cd ..
    fi
  done < $1
}

#check input.bam.list.0
# check input.bam.list.2
check input.bam.list.3
check input.bam.list.brca_wxsInwgsDownloadlist2AND3.bam.tumor_02052014.txt_part1_10
check input.bam.list.brca_wxsInwgsDownloadlist2AND3.bam.tumor_02052014.txt_part2_10

##submit job for final vcf file
#/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall /ifs/scratch/c2b2/TCGA/data/BRCA/WXS/TCGA-A1-A0SD-01A-11D-A10Y-09.bam
