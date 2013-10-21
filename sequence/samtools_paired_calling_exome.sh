#!/bin/bash
#$ -cwd

uname -a
date
## Script to perform joint calling on a list of input BAMs. Requires a file with the full paths of input files. Also requires a global settings file.

INP=""
TEMP=""
MEM="4"
AUTO="no"
OUT=""

USAGE="Usage: $0 -I <File containing list of all input BAMs> -g <global_settings.sh file> -O <path/output.vcf :optional> -t <temp dir:optional> -m <memory for job in Gb : optional> -A <no parameter : optional: use -A to trigger downstream steps automatically> -h (help)"

while getopts I:O:g:t:m:h:A o
do      case "$o" in
        I)      INP="$OPTARG";;
	O)	OUT="$OPTARG";;
        g)      GLOBAL="$OPTARG";;
        t)      TEMP="$OPTARG";;
        m)      MEM="$OPTARG";;
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
	OUT="joint.vcf"
fi
if [[ -e $OUT ]]
then
	echo "$OUT already exists.. overwriting"
fi

count=1
while IFS=$'\t' read -r -a line   
do
 
target="${line[0]}:${line[1]}-${line[3]}"
echo '#!/bin/bash'  > $out.$count.sh
echo '#$ -cwd' >> $out.$count.sh
echo 'uname -a' >> $out.$count.sh
echo " $SAMTOOLS  mpileup -DSEugd 400 -q 1 -C 50 -r $target -f $REF $INPLIST | $BCFTOOLS view -vcgT pair - | vcfutils.pl varFilter -D 2000 -1 0.0000000001 -2 0.0000000001 -3 0.0000000001  > $OUT.part-$count.vcf " >> $out.$count.sh



let count=$count+1
done < $ExonFile

exit

if [[ $AUTO == "yes" ]]
then
#	CMD=" qsub -l mem={$MEM}g,time=10:: -N GATKAnnotate $BPATH/gatk_annotate_variants.sh  -I $OUT "
	CMD=" qsub -o $OUT.log.varann.o -e $OUT.log.varann.e -l mem=5G,time=5:: $BPATH/gatk_annotate_variants.sh $OUT $INP "
	echo $CMD
	$CMD
fi

date


