#!/bin/bash
#$ -cwd

#Author:J.He
#Usage: samtools_call_exomSMs.sh  <PID-SampleCode Tumor> <PID-SampleCode Control>
#Examples: samtools_call_exomSMs.sh PAYYEP-04 PAYYEP-14
#Description: this is called by do_samtools_calvar.sh, make sure all bam files are in the same working folder as the scripts 
#input:<string:full path of bam1 > <string:full path of bam2>
#output:<file:vcf file for each line in input; log file>
##--------------------------------------
#environment variables
MYREF='/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa'
EXOBED='/ifs/scratch/c2b2/ac_lab/jh3283/ref/Exome_Targeted_Regions.BED'
BAM1=$1
BAM2=$2
wd=`pwd`

if [[ ! -d $wd/log ]]
then 
  mkdir $wd/log
fi
logd=`pwd`"/log"

echo "Reference: "$MYREF
echo "Exome Bed file: "$EXOBED
echo "BAM file1:"$BAM1
echo "BAM file2:"$BAM2
echo "log dir:" $logd
bamname1=`echo $BAM1| awk 'BEGIN{FS="/"}{split($NF,a,".");print a[1]}'`
bamname2=`echo $BAM2| awk 'BEGIN{FS="/"}{split($NF,a,".");print a[1]}'`

##-------------------------------------
#check BAM file
if [[ ! -f $1 || ! -f $2 ]] ;then  
  echo "Input BAM error!" 
  exit 1 
fi

###------------------------------call raw variants
#samtools mpileup -R -DSEugd 400 -q 1 -C 50 -f $MYREF -l $EXOBED $BAM1 $BAM2  | bcftools view -vcgT pair -p 1.1 - | vcfutils.pl fillac > $wd/${bamname1}"_"${bamname2}.var.vcf

#if [[ $? == 0 ]]
#then
#  echo "Calling $bamname1_$bamname2 DONE!" 
#  exit 0
#else
#  echo "Calling $bamname1_$bamname2 ERROR!"
#  exit 1
#fi


#gatk annotate
cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/do_annotate_gatk_annovar.sh $wd/${bamname1}"_"${bamname2}.var.vcf $BAM1 $BAM2"
echo $cmd
$cmd
#
#if [[ $? == 0 ]]
#then
#  echo "GATK ANNOTATE $bamname1_$bamname2 DONE!" 
#  exit 0
#else
#  echo "GATK ANNOTATE $bamname1_$bamname2 ERROR!"
#  exit 1
#fi

#annovar
cmd="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/do_annovar_all.sh ${wd}/${bamname1}"_"${bamname2}.var.vcf.gatk.vcf ${wd}/${bamname1}"_"${bamname2}.var.vcf.gatk.vcf.annovar.vcf"
echo $cmd
$cmd

#if [[ $? == 0 ]]
#then
#  echo "ANNOVAR $bamname1_$bamname2 DONE!" 
#  exit 0
#else
#  echo "ANNOVAR $bamname1_$bamname2 ERROR!"
#  exit 1
#fi
#

#filtering
cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/do_filtering_gatk_annovar.sh ${wd}/${bamname1}"_"${bamname2}.var.vcf.gatk.vcf.annovar.vcf -E"
echo $cmd
$cmd
#if [[ $? == 0 ]]
#then
#  echo "FILTERING $bamname1_$bamname2 DONE!" 
#  exit 0
#else
#  echo "FILTERING $bamname1_$bamname2 ERROR!"
#  exit 1
#fi
#
echo "----END-----" 
echo `date`
