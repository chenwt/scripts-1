#!/bin/bash
#By: J.He
#TODO: 



runAllGatkAnno(){
  cwd=`pwd`
  tempDir=$1
  cd $tempDir
  cnt=1
  for list in `ls split_*bam.var.vcf.list`
  do
    tempvcf=`awk '{gsub(".list","",$1);print $1}'`
    echo $tempvcf
    tempbam=`awk '{gsub("var.vcf.list","",$1);print $1}'`
    echo $tempbam
    #/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-3_doGatkAnno.sh $tempDir $tempvcf $tempbam
    ((cnt +=1 ))
  done
  echo "$cnt annotation job submitted"
}
resubSplitCall(){
    tempDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-A8-A08B-01A-11D-A19H-09.bam/
    bam=TCGA-A8-A08B-01A-11D-A19H-09.bam
    cntline=0
    for cnt in  9 96 97 98 99 
    do
      ((cntline+=1))
      echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/samtools_callSomaticVar.sh $tempDir/split_$cnt.$bam" |qsub -l mem=4g,time=8:: -N reCall$cnt -e $tempDir/log -o $tempDir/log -cwd
    done
    echo $cntline call job resubmitted
}



##------------executing------
resubSplitCall 
# runAllGatkAnno /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-A8-A08B-01A-11D-A19H-09.bam/
