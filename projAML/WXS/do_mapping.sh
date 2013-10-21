#!/bin/bash
#$ -cwd
uname -a
echo "Start $0:`date`"

fq=$1	#Takes single filename not a lit	
setting=$2
automated=$3	#if initiated by automated pipeline then input argument must be ProjectID, this triggers automated downstream steps

BPATH="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome"

if [[ $setting == "" ]]; then
    setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Exome/global_setting_b37.sh"
fi

. $setting 

#for f in `cat $fq1list`  ## list.forwardReads.fq.txt`
#  do 
# example: R_U_NI_D_laneX_SAMPLEID_INDEX_1.fastq
  job_id=`basename $fq |  sed 's/\_1.fastq//' `	#unique name for jobname (includes runid_lane#_sampleID
  g=`echo $fq | sed 's/_1.fastq/_3.fastq/'`
  nthread=2
  if  [ ! -e $g ];  then g="" ;  nthread=2 ;else   nthread=4;fi

sampleName=`basename $fq | sed 's/\_1.fastq//' | cut -f6  -d '_'`	#acc to  new convention, to change to old do -f4
rg1=`basename $fq|cut -f3 -d"_"`	#the 3rd field in the RUNID
ln=`basename $fq|cut -f5 -d"_" | sed 's/lane//'`	#lane#
readgroup=$rg1"_"$ln"_"$sampleName	#unique Readgroup based on library
ID="$rg1$ln"				#shorter ID based on readgroup

if [ ! -d mapping"$ln" ]; then mkdir mapping"$ln" ; fi
if [ ! -d logs"$ln" ]; then mkdir logs"$ln" ; fi
  
if [[ $automated == "" ]]; #was NOT triggered by automatic pipeline
then
  cmd="qsub -pe smp $nthread -R y -l mem=4G,time=48:: -o logs$ln/mapping.$sampleName.o -e logs$ln/mapping.$sampleName.e -N map$ln.$job_id $BPATH/mapping-two-cores.sh -i $fq -y $ID -z $readgroup -n $sampleName -s $setting -o mapping$ln/$sampleName -t 2 -c 1 -p $g  "
else	#trigger automatic process
  cmd="qsub -pe smp $nthread -R y -l mem=4G,time=48:: -o logs$ln/mapping.$sampleName.o -e logs$ln/mapping.$sampleName.e -N map$ln.$job_id.AUTO $BPATH/mapping-two-cores.sh -i $fq -y $ID -z $readgroup -n $sampleName -s $setting -o mapping$ln/$sampleName -t 2 -c 1 -A $automated -p $g "
fi
  echo $cmd >> logs$ln/history.$sampleName.txt
  echo "Qsub :`date`" >> logs$ln/history.$sampleName.txt
  $cmd
#done

echo "End $0:`date`"

