#!/bin/bash
#By: J.He
#TODO: more parameter checking 
##Desp.: invoked in step-1_splitBAM.sh and invoke step-3_gatkAnno.sh
##input: splitted bam file folder; subBam full path

##---setting parameters-
sh /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh 
srcDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
tempDir=$1
subBAM=$2
ERR="[ERR:]"
MSG="[MSG:]"
DONE="[DONE:]"
rootd=`pwd`

subBamName=`echo $subBAM|awk 'BEGIN{FS="/"}{print $NF}'`
cnt=`echo $subBamName|awk -F. '{gsub("split_","",$1);print $1} '`
pid=`echo $subBamName|awk 'BEGIN{FS="."}{split($2,a,"-");print a[3]} '`

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

echo "$srcDIR/samtools_callSomaticVar.sh $subBAM" |qsub -N c_${pid}_$cnt -e $logDir -o $logDir -l mem=8G,time=40:: -cwd >>$cwd/qsubCall.log
tail -1 $cwd/qsubCall.log

##-----invoke annoatation
cd $cwd
if [ ! -d $cwd/annoJobSubmitted ]; then mkdir $cwd/annoJobSubmitted ; fi 
cntNF=`ls -1 $subBAM.var.vcf.list 2>/dev/null | wc -l`
cntjob=0
echo -e "new var.vcf file\t" $cntNF
while [ $cntNF -le 1 ] 
do
    if [ $cntNF -eq 1 ] ; then
         subVCF=`readlink -f $subBAM.var.vcf.list|awk '{gsub(".list","",$1);print $1}'`
         cmd="$srcDIR/step-3_doGatkAnno.sh $cwd $subVCF $subBAM"	
	 $cmd 
	 subVCFname=`echo $subVCF|awk -F"/" '{print $NF}'`
	 mv $subBAM.var.vcf.list $cwd/annoJobSubmitted/
	 break
    fi
    echo "waiting calling file...$(date)"
    cntNF=`ls -1 split_*bam.var.vcf.list 2>/dev/null | wc -l`
    sleep 10m
done

echo "$subBAM anotation job submitted!"
echo "#---DONE----start-annotation-----"
