#!/bin/bash
#$ -cwd
# Modified by J.HE: Oct 7th
# Input: <samtools output vcf file>, <original bam file> 
# Output: <annotated VCF file> 

uname -a
msg="<MESSAGE>: "
cmd="<COMMAND>: "
warn="<WARNING>: "
err="<ERROR>: "	
echo "$msg Start $0: `date`"

USAGE="--USAGE : $0 <INPUT.vcf> <INPUT.bam> "
VCF=$1
BAM=$2
GLOBAL=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh
. $GLOBAL

if [[ $VCF == "" || $BAM == "" || ! -s $VCF ]]
then
        echo $USAGE
        exit 1
fi

vcfName=$(echo $VCF|awk 'BEGIN{FS="/"}{print $NF}')
bamName=$(echo $BAM|awk 'BEGIN{FS="/"}{print $NF}')


OUTPUT=$VCF".gatk.vcf" 
TEMP=$OUTPUT"_temp"

if [ ! -d $TEMP ]; then mkdir -p $TEMP ; fi

echo $SAMTOOLS
if [ ! -f $BAM.bai ];then 
  echo "index bam...."
  $SAMTOOLS index $BAM 
fi
echo "annotating..." 	
## List of annotations from /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/Exome/Exome_newGATK/Variant_Annotation_List.txt.new 
java -Xmx3g -Djava.io.tmpdir=$TEMP -jar $GATKJAR \
	-R $REF \
	-T VariantAnnotator  \
	--variant $VCF \
	-o $OUTPUT \
	-I $BAM \
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
	rm -R $TEMP
	echo "$msg End $0: `date`"
	exit 0
else
	echo "$msg End $0: `date`"
	exit 1
fi

echo "GATKannotationDONE" > $VCF.gatk.vcf.list
