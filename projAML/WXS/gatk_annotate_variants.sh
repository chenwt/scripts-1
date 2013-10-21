#!/bin/bash
#$ -cwd
# Modified by J.HE: Oct 7th
# Input:
# Output:

uname -a
msg="<MESSAGE>: "
cmd="<COMMAND>: "
warn="<WARNING>: "
err="<ERROR>: "	
echo "$msg Start $0: `date`"
	
exome="NO"
VCF=$1
BAM=$2
GLOBAL=$3


# while getopts v:b:s:c:h:eA o
# do      case "$o" in
#        	v)      VCF="$OPTARG";;
# 		b)	BAM="$OPTARG";;
# 		e)	exome="YES";;
# 		c)	CHR="$OPTARG";;
# 		s)  GLOBAL="$OPTARG";;
# 		A)	AUTO="yes";;
#       	h)  echo $USAGE
#                exit 1;;
#        esac
# done
# #
echo $VCF
echo $GLOBAL
echo $BAM
BAM1=$(awk '{print $1}' $BAM)
BAM2=$(awk '{print $2}' $BAM)
# bam list contains two files

USAGE="--USAGE : $0 <INPUT.vcf> <INPUT.bam.list> <YES/NO> <GLOBAL.setting.file> "

if [[ $VCF == "" || $GLOBAL == "" || $BAM == "" || ! -s $VCF ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL


## VCF=$1		# VCF="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/SAMPLES/AC2/VarCalling_samtool/1.var.flt.vcf"
OUTPUT=$VCF".gatk.vcf" 
## BAM=$2  	##Argument 2 is file with list of bams.

if [[ $exome == "YES" ]]
then
 	REF="/ifs/data/c2b2/ngs_lab/ngs/usr/src/bwa-0.7.3/human_b37/human_g1k_v37.fasta"
	## i.e. this script was called from nextera exome pipeline
fi
 
TEMP=$OUTPUT"_temp"
if [ ! -e $TEMP ]; then mkdir $TEMP ; fi

if [ ! -f $BAM1.bai ];then $SAMTOOLS index $BAM1 ; fi
if [ ! -f $BAM2.bai ];then $SAMTOOLS index $BAM2 ; fi
	
# List of annotations from /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/Exome/Exome_newGATK/Variant_Annotation_List.txt.new 
java -Xmx3g -Djava.io.tmpdir=$TEMP -jar $GATKJAR \
	-R $REF \
	-T VariantAnnotator  \
	--variant $VCF \
	-o $OUTPUT \
	-I $BAM1 \
	-I $BAM2 \
	-L $VCF \
	--dbsnp $DBSNP132 \
	-XA "SnpEff" \
	-XA "TransmissionDisequilibriumTest" \
	-XA "ChromosomeCounts" \
	-XA "HardyWeinberg" \
	-XA "NBaseCount" \
	-XA "BaseCounts" \
	-XA "MVLikelihoodRatio" \
	-XA "RodRequiringAnnotation" \
	-XA "TechnologyComposition" \
	-XA "SampleList" \
	-all \
    -A DepthPerAlleleBySample  -A AlleleBalance  -A DepthOfCoverage  -A BaseQualityRankSumTest  -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth  -A RMSMappingQuality  -A SpanningDeletions  -A HaplotypeScore  -rf BadCigar 
## -filterMBQ


if [[ $? == 0 ]]
then
	echo "$msg End $0: `date`"
	exit 0
else
	echo "$msg End $0: `date`"
	exit 1
fi

