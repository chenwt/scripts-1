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
	A)	AUTO="no";;
        h)      echo $USAGE
                exit 1;;
        esac
done
echo "INPUT:" $INP
echo "GLOBAL:"$GLOBAL
echo "OUT:" $OUT
echo "AUTO:" $AUTO

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

CMD=" $SAMTOOLS  mpileup -DSEugd 400 -q 1 -C 50 -f $REF $INPLIST | $BCFTOOLS view -p 1.1 -vcg - | vcfutils.pl fillac - >$OUT" 
# CMD=" $SAMTOOLS  mpileup -DSEugd 400 -q 1 -C 50 -f $REF $INPLIST | $BCFTOOLS view -p 1.1 -vcg - | vcfutils.pl varFilter -D 2000 -1 0.0000000001 -2 0.0000000001 -3 0.0000000001 " 
echo $CMD
echo "output: $OUT "
 $SAMTOOLS  mpileup -DSEugd 400 -q 1 -C 50 -f $REF $INPLIST | $BCFTOOLS view -p 1.1 -vcg - | vcfutils.pl fillac - >$OUT

# if [[ $AUTO == "yes" ]]
# then
#         cd ${INP}
#         echo $(pwd)

#         if [ ! -d log ]; then
#                 echo " making log dir"
#                 mkdir log
#         fi
#         logd=$(readlink -m ./log)
#         echo ${logd}
#         for chr in $(seq 1 22) X Y
#         do 
#                 cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/gatk_annotate_variants.sh \
#                 ${INP}/out_pair${chr}.vcf  \
#                 ${INP}/bams${chr}.list \
#                 /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/global_setting_GRC37lite.sh"
#                 # echo ${cmd}
#                 echo ${cmd} | qsub -l mem=4g,time=10:: -N annot${chr} -e ${logd}/ -o ${logd}/  >> ${logd}/qsub_log
#                 tail -1  ${logd}/qsub_log
#         done 


# #	CMD=" qsub -l mem={$MEM}g,time=10:: -N GATKAnnotate $BPATH/gatk_annotate_variants.sh  -I $OUT "
# 	# CMD=" qsub -o $OUT.log.gatk.o -e $OUT.log.gatk.e -l mem=5G,time=5:: $BPATH/gatk_annotate_variants.sh $OUT $INP $GLOBAL "
# 	# echo $CMD
# 	# $CMD
# fi

date
