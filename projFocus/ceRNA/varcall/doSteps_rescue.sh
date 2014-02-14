#!/bin/bash
#By: J.He
#$-cwd
#TODO: 


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

runAllCall(){
  cwd=`pwd`
  tempDir=$1
  cd $tempDir
  cnt=1
  # for list in `ls split_*bam.list`
  for tempbam in `ls split_*bam`
  do
    # tempbam=`echo $list|awk '{gsub(".list","",$1);print $1}'`
    # vcf=$tempbam.var.vcf.list
     # echo $tempbam
     if [ -f $tempbam ]; then
	/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh $tempDir $tempbam
	((cnt +=1 ))
     fi 
  done
  echo "$cnt calling job submitted"
}

runAllGatkAnno(){
  cwd=`pwd`
  tempDir=$1
  cd $tempDir
  cnt=1
  for list in `ls split_*bam.var.vcf`
  do
    # tempvcf=`echo $list|awk '{gsub(".list","",$1);print $1}'`
    tempvcf=$list
     # echo $tempvcf
    tempbam=`echo $list|awk '{gsub(".var.vcf","",$1);print $1}'`
     # echo $tempbam
     if [[ -f $tempvcf ]] && [[ -f $tempbam ]] ; then 
       /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-3_doGatkAnno.sh $tempDir $tempvcf $tempbam
       ((cnt +=1 ))
     fi
  done
  echo "$cnt annotation job submitted"
}

##------------executing------
## resubSplitCall 
# runAllGatkAnno /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-A8-A08B-01A-11D-A19H-09.bam/
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-A8-A08B-01A-11D-A19H-09.bam /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/  
#---one sample try out
# qsub -l mem=6g,time=24:: -N t_A0J6_103  -cwd  /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/temp-region_103.sh
# qsub -l mem=8g,time=24:: -N t_A0J6_103  -cwd  /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/temp-region_102.sh
# qsub -l mem=8g,time=24:: -e log/ -o log/ -N ta_A0J6_100  -cwd  /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/temp-region_100.sh
# for cnt in `seq 1 44`
# # for cnt in `seq 45 50`
# # for cnt in `seq 50 102`
# do
# qsub -l mem=8g,time=24:: -e log/ -o log/ -N t_A0J6_$cnt -cwd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/temp-region_$cnt.sh
# done 

# for cnt in `seq 1 50` 101 102 103
# do
#   rm /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/split_$cnt.TCGA-AO-A0J6-01A-11D-A128-09.bam.var.vcf.list
# done
# runAllCall /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/
#/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/split_103.TCGA-AO-A0J6-01A-11D-A128-09.bam

# runAllGatkAnno /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam 
# runAllGatkAnno /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-BH-A0BM-01A-11D-A060-02.bam  

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/split_101.TCGA-AO-A0J6-01A-11D-A128-09.bam

# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam/split_102.TCGA-AO-A0J6-01A-11D-A128-09.bam

# tempDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-BH-A0BM-01A-11D-A060-02.bam
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/  

# tempDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-AO-A0J6-01A-11D-A128-09.bam
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-3_doGatkAnno.sh  $tempDir $tempDir/split_102.TCGA-AO-A0J6-01A-11D-A128-09.bam.var.vcf $tempDir/split_102.TCGA-AO-A0J6-01A-11D-A128-09.bam
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-3_doGatkAnno.sh  $tempDir $tempDir/split_101.TCGA-AO-A0J6-01A-11D-A128-09.bam.var.vcf $tempDir/split_101.TCGA-AO-A0J6-01A-11D-A128-09.bam
 # /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/  

runAllCall /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-A8-A08L-01A-11D-A19H-09.bam 
