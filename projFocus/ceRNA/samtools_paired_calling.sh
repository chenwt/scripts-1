#!/bin/bash
#$ -cwd

uname -a
echo "Start $0: `date` "
## Script to perform paired calling on BAM pairs. 
# Requires:     1. <a file with the full paths of input files>
#               2. <a global settings file>
# ~/scripts/projFocus/ceRNA/samtools_paired_calling.sh  \\ 
#  -I bam.list \\
#  -g ~/scripts/projFocus/ceRNA/global_setting_GRC37lite.sh \\
#  -O /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/result.vcf \\
# optional below
#  -t <temp dir:optional> 
#  -m 4
#  -A no

INP=""
MEM="8"
AUTO="no"
OUT=""

USAGE="Usage: $0 -I <File containing list of all input BAMs> -g <global_settings.sh file> -O <path/output.vcf :optional> -t <temp dir:optional> -m <memory for job in Gb : optional> -A <no parameter : optional: use -A to trigger downstream steps automatically> -h (help)"

while getopts I:O:g:m:h:A o
do      case "$o" in
        I)      INP="$OPTARG";;
	O)	OUT="$OPTARG";;
        g)      GLOBAL="$OPTARG";;
        m)      MEM="$OPTARG";;
	A)	AUTO="yes";;
        h)      echo $USAGE
                exit 1;;
        esac
done

# echo "INPUT:" $INP
# echo "GLOBAL:"$GLOBAL
# echo "OUT:" $OUT
# echo "AUTO:" $AUTO


if [[ $INP == "" || $GLOBAL == "" || ! -s $INP ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

MYSCDIR="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/"
# echo $MYSCDIR

INPLIST=`cat $INP`
if [[ $OUT == "" ]]
then
        OUT="joint.vcf"
fi
if [[ -e $OUT ]]
then
        echo "$OUT already exists.. overwriting"
fi

CMD="${SAMTOOLS} mpileup -DSEugd 400 -q 1 -C 50 -f $REF $INPLIST | $BCFTOOLS view -vcgT pair -p 1.1 - | vcfutils.pl fillac -  > $OUT "
echo $CMD 
# echo "output: $OUT " 
$SAMTOOLS mpileup -DSEugd 400 -q 1 -C 50 -f $REF $INPLIST | $BCFTOOLS view -vcgT pair -p 1.1  - | vcfutils.pl fillac -  > $OUT
# #try to output allele frequency
# #$SAMTOOLS mpileup -DSEugd 400 -q 1 -C 50 -f $REF $INPLIST | $BCFTOOLS view -vcgT pair - | vcfutils.pl varFilter -D 2000 -1 0.0000000001 -2 0.0000000001 -3 0.0000000001  > $OUT

##to runn the following steps, first index all splitted bam files
if [[ $AUTO == "yes" ]]
then
        INPD=$(echo ${INP}|awk 'BEGIN{OFS=FS="/"}{$NF="";print $0}')
        cd ${INPD}
        echo $(pwd)

        if [ ! -d log ]; then
                echo " making log dir"
                mkdir log
        fi
        logd=$(readlink -m ./log)
        # echo ${logd}
        for chr in $(seq 1 22) X Y
        do 
                cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/gatk_annotate_variants.sh \
                ${INPD}/out_pair${chr}.vcf  \
                ${INPD}/bams${chr}.list \
                /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/global_setting_GRC37lite.sh"
                # echo ${cmd}
                echo ${cmd} | qsub -l mem=4g,time=10:: -N annot${chr} -e ${logd}/ -o ${logd}/  >> ${logd}/qsub_log
                tail -1  ${logd}/qsub_log
        done 

# ##	CMD=" qsub -l mem={$MEM}g,time=10:: -N GATKAnnotate $BPATH/gatk_annotate_variants.sh  -I $OUT "
#  #       # CMD=" qsub -o $OUT.log.varann.o -e $OUT.log.varann.e -l mem=5G,time=5:: $BPATH/gatk_annotate_variants.sh $OUT $INP"
#  ##        CMD=" qsub -o $OUT.log.gatk.o -e $OUT.log.gatk.e -l mem=8G,time=5:: ${MYSCDIR}/gatk_annotate_variants.sh $OUT $INP $GLOBAL"
#  # #       echo $CMD
# 	# # CMD=" qsub -o $OUT.log.anno.o -e $OUT.log.anno.e -l mem=5G,time=5:: sh $ANNOVAR/do_annovar_all.sh $OUT annovar_$OUT"
# #	# $CMD

fi

echo "End $0 `date` "

exit

