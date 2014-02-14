#!/bin/bash
#$ -cwd

uname -a
msg="<MESSAGE>: "
cmd="<COMMAND>: "
warn="<WARNING>: "
err="<ERROR>: "	
echo "$msg Start $0: `date`"
## Script to perform single calling  on a list of input BAMs. Requires a file with the full paths of input files. Also requires a global settings file.

INP=""
AUTO="no"
OUT=""
CHR=""

USAGE="Usage: $0 -I <File containing list of all input BAMs> -g <global_settings.sh file> -O <path/output.vcf :optional> -A <flag to trigger AUTO annotation> -h <flag> HELP"

while getopts I:O:N:D:C:g:hA o
do      case "$o" in
        I)      INP="$OPTARG";;
	O)	OUT="$OPTARG";;
	N)	sampleName="$OPTARG";;
	T)	Tumor="$OPTARG";;
	D)	Dir="$OPTARG";;
	C)	CHR="$OPTARG";;
	g)      GLOBAL="$OPTARG";;
	A)	AUTO="yes";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $INP == "" || $GLOBAL == "" || ! -s $INP ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

INPLIST=`cat $INP`
if [[ $OUT == "" ]]
then
	OUT="single.vcf"
fi
if [[ -e $OUT ]]
then
	echo "$OUT already exists.. overwriting"
fi

if [[ $CHR == "" ]]
then 
	ropt=""
else
	ropt=" -r $CHR "
fi

CMD=" $SAMTOOLS  mpileup -DSEugd 400 -q 1 -C 50 $ropt -f $REF $INPLIST | $BCFTOOLS view -p 0.9 -vcg - "
## | vcfutils.pl varFilter -D 2000 -1 0.0000000001 -2 0.0000000001 -3 0.0000000001 "
echo $CMD 
echo "output: $OUT " 
$SAMTOOLS  mpileup -DSEugd 400 -q 1 -C 50 $ropt -f $REF $INPLIST | $BCFTOOLS view -p 0.9 -vcg - > $OUT 
## | vcfutils.pl varFilter -D 2000 -1 0.0000000001 -2 0.0000000001 -3 0.0000000001  > $OUT

if [[ $? == 0 ]]
then
	echo "Done Calling"
	if [[ $AUTO == "yes" ]]
	then
#	CMD=" qsub -o $OUT.log.varann.o -e $OUT.log.varann.e -l mem=5G,time=5:: $BPATH/gatk_annotate_variants.sh -v $OUT -b $INP -s $GLOBAL -c $CHR"
		CMD="sh $BPATH/gatk_annotate_variants.sh -v $OUT -b $INP -s $GLOBAL -c $CHR "
		echo $CMD
		$CMD 1> $OUT.log.varann.o  2> OUT.log.varann.e 
		if [[ $? == 0 ]]
		then
			echo "Success GATK-ANNOTATE"
			#Do ANNOVAR genomic annotation
			CMD="sh $BPATH/do_annovar_all.sh $OUT "
			echo $CMD
			$CMD 1> $OUT.log.annovar.o  2> $OUT.log.annovar.e
			if [[ $? == 0 ]]
			then
				echo "Success ANNOVAR"
				echo "End $0: `date`"
				exit 0
			else
				echo "Failed ANNOVAR"
				echo "End $0: `date`"
				exit 1
			fi
		else
			echo "Failed GATK-ANNOTATE"
			echo "End $0: `date`"
			exit 1
		fi
	else
		echo "$msg End $0: `date`"
		exit 0
	fi
else
	echo "Failed Calling"
	echo "End $0 `date` "
	exit 1
fi


