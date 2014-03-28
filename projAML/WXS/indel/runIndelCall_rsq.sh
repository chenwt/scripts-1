#!/bin/bash
#By: J.He
#TODO: 

#data folder:
dataDir=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/
outDir=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/indel/

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

###debugging for first run
###waiting for all windown file job finished 
check() {
  while read line
  do
    bam=`echo $line|awk 'BEGIN{FS="/"}{print $NF}'` 
    if [ -d temp-$bam ] ;then
      cd temp-$bam
      jobFinished=`awk 'END{print NR}' jobFinished.log `
      jobSubmitted=`ls win*sh 2>/dev/null|wc -l`
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

# check done/pid.txt.inputBamlist_1
# check pid.txt.inputBamlist_2
# check pid.txt.inputBamlist_3
# check pid.txt.inputBamlist_4

##submit job for final vcf file

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PASFEW/tumor_sample/TARGET-20-PASFEW-09A-01D/SRX100073/SRR564140/SRR564140.rmdup.new.bam  
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PANSBH/normal_sample/TARGET-20-PANSBH-14A-01D/SRX100041/SRR563909/SRR563909.rmdup.new.bam 

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PARGVC/normal_sample/TARGET-20-PARGVC-14A-01D/SRX100062/SRR563748/SRR563748.rmdup.new.bam

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PARGVC/tumor_sample/TARGET-20-PARGVC-04A-01D/SRX100061/SRR563639/SRR563639.rmdup.new.bam 

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PANSBH/normal_sample/TARGET-20-PANSBH-14A-01D/SRX100041/SRR563909/SRR563909.rmdup.new.bam 

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PANTWV/tumor_sample/TARGET-20-PANTWV-04A-01D/SRX100043/SRR563879/SRR563879.rmdup.new.bam 

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir /ifs/scratch/c2b2/TCGA/data/AML/29212/reads/AML/TARGET-20-PANTWV/normal_sample/TARGET-20-PANTWV-14A-01D/SRX100010/SRR563910/SRR563910.rmdup.new.bam 

for line in `ls -1 |grep temp|awk -F- '{print $2}'`
do
    bamPath=`grep $line pid.txt.inputBamlist`
    # echo $line $bamPath
    /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3reqsub.sh $outDir $bamPath 
done


