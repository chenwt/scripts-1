#!/bin/bash
#By: J.He
#TODO: 
##Desp.:singel call variants
##input: splitted bam file folder; subBam full path

##---setting parameters-
sh /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh 
# REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/humanGenome_hg19_regions104.bed
REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/test/hg19_region_test.bed
srcDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
tempDir=$1
subBAM=$2
ERR="[ERR:]"
MSG="[MSG:]"
DONE="[DONE:]"
rootd=`pwd`
###-----test----
# subBAM=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/test/temp-TCGA-A8-A08B-01A-11D-A19H-09.bam/split_1.TCGA-A8-A08B-01A-11D-A19H-09.bam
##-----test--end
# echo $subBAM
subBamName=`echo $subBAM|awk 'BEGIN{FS="/"}{print $NF}'`
# echo $subBamName
cnt=`echo $subBamName|awk -F. '{gsub("split_","",$1);print $1} '`
pid=`echo $subBamName|awk 'BEGIN{FS="."}{split($2,a,"-");print a[3]} '`
# echo $cnt
# echo $pid
###-----check dir setting-----
if [ ! -d $tempDir ]; then 
  echo -e $ERR"temp dir not there!" 
  exit
fi
  
if [ ! -d $tempDir/log ]; then 
  echo -e $ERR"log dir not there!" 
  exit
fi

cwd=$tempDir
logDir=$cwd/log

cd $cwd
###----call vars
# echo -e "calling on \t $subBAM"   
# echo -e "current working dir\t $cwd"

# cp $srcDIR/samtools_callSomaticVar.sh .
# echo "$srcDIR/samtools_callSomaticVar.sh $subBAM" |qsub -N c_${pid}_$cnt -e $logDir -o $logDir -l mem=8G,time=40:: -cwd >>$cwd/qsubCall.log
# tail -1 $cwd/qsubCall.log


##-----invoke annoatation
cd $cwd
if [ ! -d $cwd/annoJobSubmitted ]; then mkdir $cwd/annoJobSubmitted ; fi 
sleep 2s
# sleep 30m
cntNF=`ls -1 $subBAM.var.vcf.list 2>/dev/null | wc -l`
cntjob=0
echo "new file.." $cntNF
while [ $cntNF -le 1 ] 
do
    if [ $cntNF -eq 1 ] ; then
         subVCF=`readlink -f $subBAM.var.vcf.list|awk '{gsub(".list","",$1);print $1}'`
         cmd="$srcDIR/test_step-3_doGatkAnno.sh $cwd $subVCF $subBAM"	
	 $cmd 
	 subVCFname=`echo $subVCF|awk -F"/" '{print $NF}'`
         # echo "$subVCF annotating started" >> $cwd/annoJobSubmitted/anno_${subVCFname} 
	 mv $subBAM.var.vcf.list $cwd/annoJobSubmitted/
	 break
    fi
    echo "waiting calling file...$(date)"
    cntNF=`ls -1 split_*bam.var.vcf.list 2>/dev/null | wc -l`
    sleep 1s
    # sleep 10m
done

echo "$subBAM anotation job submitted!"
echo "#---DONE----start-annotation-----"
