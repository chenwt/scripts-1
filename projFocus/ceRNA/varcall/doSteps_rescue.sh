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
    tempbam=`echo $list|awk '{gsub(".var.vcf","",$1);print $1}'`
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

# runAllCall /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/temp-TCGA-A8-A08L-01A-11D-A19H-09.bam 

 checkSubBAM(){
   tempDir=$1
   echo $tempDir
   BAM=`echo $tempDir|awk -F"/" '{gsub("temp-","",$NF);print $NF}'`
   ls $tempDir/split_*.$BAM.list 2>/dev/null |wc -l
   ls $tempDir/split_*.$BAM.var.vcf.list 2>/dev/null |wc -l
   ls $tempDir/split_*.$BAM.var.vcf.gatk.vcf.list 2>/dev/null |wc -l
 }

 function doRunGATKAnnoAll(){
  ###-----------do annotation ------
  
  tempDir_array=( `ls -r |grep "temp" ` )
  # echo ${tempDir_array[@]}
  for tempDir in "${tempDir_array[@]}"
  do
      # echo $tempDir
      cntVCF=`ls $tempDir/*var.vcf 2>/dev/null |wc -l`
      cntBAM=`ls $tempDir/split*.bam  2>/dev/null |wc -l`
      # cntTotal=`wc -l $REGIONFILE 2>/dev/null|wc -l `
      cntTotal=103
      echo "vcf $cntVCF bam $cntBAM total $cntTotal" 
      if [ $cntBAM != $cntTotal ] || [ $cntVCF != $cntTotal ]
      then
          echo -e "skip Annotation $tempDir "
	  continue 
      fi
      echo "Annotation $tempDir" 
      pid=`echo $tempDir|awk -F"-" '{print $4}' `
      dir=`readlink -m $tempDir`
      runAllGatkAnno $dir  
  done 
}
runReQsub(){
  tempDir=$1
  for tempBAM in `ls $tempDir/split*bam`
  do
    cnt=`echo $tempBAM|awk -F"/" '{split($NF,a,".");gsub("split_","",a[1]);print a[1]}'`
    pid=`echo $tempDir|awk -F"/" '{split($NF,a,"-");print a[4]}'`
    echo $cnt" " $pid 
    qsub -l mem=4g,time=48:: -N ${pid}_${cnt} -cwd -e $tempDir/log -o $tempDir/log $tempDir/temp-region_$cnt.sh 
  done
}
rootDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs

# runReQsub $rootDir/temp-TCGA-C8-A130-01A-31D-A19H-09.bam
# runReQsub $rootDir/temp-TCGA-A2-A04Q-01A-21D-A128-09.bam
# runReQsub $rootDir/temp-TCGA-AR-A1AY-01A-21D-A12Q-09.bam
# runReQsub $rootDir/temp-TCGA-A2-A259-01A-11D-A314-09.bam 
# runReQsub $rootDir/temp-TCGA-AN-A04D-01A-21D-A314-09.bam

# checkSubBAM  TCGA-A2-A04Q-01A-21D-A128-09.bam
# checkSubBAM  TCGA-E9-A1NH-01A-11D-A14G-09.bam

# runAllGatkAnno  $rootDir/temp-TCGA-AO-A0J4-01A-11D-A128-09.bam
# for cnt in 8 22 77 78 100 102 103
# do
   # /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh  $rootDir/temp-TCGA-A8-A08L-01A-11D-A19H-09.bam $rootDir/temp-TCGA-A8-A08L-01A-11D-A19H-09.bam/split_$cnt.TCGA-A8-A08L-01A-11D-A19H-09.bam
# done

# cnt=11
# tempBAM=TCGA-C8-A130-01A-31D-A19H-09.bam
# /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh  $rootDir/temp-$tempBAM $rootDir/temp-$tempBAM/split_$cnt.$tempBAM

# tempDir=$rootDir/temp-TCGA-AO-A0J4-01A-11D-A128-09.bam
 # /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/  
# tempDir=$rootDir/temp-TCGA-A8-A08L-01A-11D-A19H-09.bam
  # /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/  

# for cnt in `seq 70 103` 
# for cnt in `seq 1 103` 
# do
#     tempDir=$rootDir/temp-TCGA-A2-A04Q-01A-21D-A128-09.bam
#     bam=$tempDir/split_$cnt.TCGA-A2-A04Q-01A-21D-A128-09.bam
#     if [ ! -f $bam.var.vcf ] ; then
# 	/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-2_callVar.sh  $tempDir $bam
#     else
#       rm $bam
#     fi
# done

 # tempDir=$rootDir/temp-TCGA-A2-A04Q-01A-21D-A128-09.bam
 # /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/step-4_mergeVCF.sh $tempDir /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/  
# doRunGATKAnnoAll

 for temp in `ls -r |grep temp`
 do 
  checkSubBAM $temp
 done
