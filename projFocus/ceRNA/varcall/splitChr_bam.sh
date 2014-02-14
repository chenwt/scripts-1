#!/bin/bash
#$ -cwd
#$ -l mem=4G,time=4::
#Desp.: Takes a BAM file and splits it equally chunks, and also extracts all unpaired read pairs, don't process them. Output directory is split.N.<input>.
##input: <full path of input bam file>
##output: <splited bam files> <files containing the full path of splitted bam files>
#J.HE: 

USAGE="Usage: $0 -n name-sample -i /path/input.bam  -s /path/global_setting -o path/output-directory/  -h [flag to show help]"
uname -a
echo "Start $0:`date`"
starttime=`date +%s`
function doexit(){	##first argument is exit status
	endtime=`date +%s`
	runtime=$((endtime-starttime))
	exit $1
}

while getopts i:o:l:h opt
  do
  case "$opt" in
	i) input="$OPTARG";;
    o) output="$OPTARG";;
	l) LOGDIR="$OPTARG";;
    h)    echo $USAGE
          doexit 1;;
  esac
done
##-------set and check parameter
SAMTOOLS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools
READSNUMBER=15000000
sampleName=`echo $input| awk 'BEGIN{FS="/"}{print $NF}'`


echo -e "input:\t"$input
echo -e "bam file:\t"$sampleName
if [[ $input == "" $sampleName == "" ]]; then
    echo $USAGE
    doexit 1
fi

##-----create temp folder-----
input=`readlink -f $inp`
if [[ -d temp-$sampleName ]]; then mkdir temp-$sampleName; fi
if [[ -d log ]]; then mkdir log; fi
cd temp-$sampleName
cwd=`pwd`
if [[ -d log ]]; then mkdir log; fi
$LOGDIR=`pwd`/log
if [[ $output == "" ]]
then
    OUTDir=$cwd/
else
  output=`readlink -f $output `
fi
if [ ! -d $output ]; then mkdir $output ; fi

##-------Step0 Get header-----

if [[ ! -s $output/"head.sam" ]]
then
  $SAMTOOLS view -H $input > $output/"head.sam"
fi

##Step1 Extract unmapped reads also from the BAM (for geneFusion analysis)
cmdunmapped="qsub -o $LOGDIR/split.$sampleName.unmap.o -e $LOGDIR/split.$sampleName.unmap.e -l mem=2G,time=4:: -N split.$sampleName.unmapped $BPATH/get_unmapped_reads.sh -i $input -o $output/unmapped.pairs.bam -s $setting -l $LOGDIR "
echo $cmdunmapped >> $LOGDIR/split."$sampleName"."unmap".bam.log 
$cmdunmapped 

##Step 2: Split the bam by given reads number.
#example of head.sam
#@HD     VN:1.0  GO:none SO:coordinate
#@SQ     SN:1    LN:249250621
# for chr in `seq 1 22` X Y MT GL
# do
	echo "Extracting $chr"
	$SAMTOOLS view -h $sampleName |awk -v n=$READSNUMBER -v f=$sampleName 'NR%n==1{x="split"++i"_"f;}{print > x}'  
	# echo "Done split $input $chr" >> $LOGDIR/split."$sampleName"."$chr".bam.o 
	# sh  $BPATH/do_mergeChr_bam.sh -n $sampleName -s $setting -l $LOGDIR -c $chr -p 26
# done 

##----output full path the splitted files
echo -n "" > $sampleName.splitted.files.list 
for file in `ls cwd/split*$samplename`
do
  readlink -f $file >> $sampleName.splitted.files.list
done

echo "End $0: `date`"
doexit 0


