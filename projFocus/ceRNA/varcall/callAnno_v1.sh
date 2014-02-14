#!/bin/bash
#By: J.He
#$ -cwd
##------------parameters
BAM=$1

##startcopy
##--------------------
uname -a
msg="<MESSAGE>: "
cmd="<COMMAND>: "
warn="<WARNING>: "
err="<ERROR>: "	
echo "$msg calling Start $0: `date`"
USAGE="--USAGE : $0 <INPUT.bam> "

outputDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars/
cwd=`pwd`

sh /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh
srcDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
SAMTOOLS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools
BCFTOOLS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bcftools
VCFUTILS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/vcfutils.pl

##------------call
##-D output per-sample read depth -S per-sample Phre-scaled strand bias p-value -u for piping -f faidx index reference file 
$SAMTOOLS mpileup -DSuf $REF $BAM |$BCFTOOLS view -bvcg -p 0.9 - > $BAM.var.bcf
$BCFTOOLS view $BAM.var.bcf | $VCFUTILS varFilter -D 400 > $BAM.var.vcf
rm $BAM.var.bcf
echo $BAM.var.vcf > $BAM.var.vcf.list

###-----------anno
bamName=$(echo $BAM|awk 'BEGIN{FS="/"}{print $NF}')
VCF=$BAM.var.vcf
vcfName=$(echo $VCF|awk 'BEGIN{FS="/"}{print $NF}')
OUTPUT=$VCF".gatk.vcf" 
TEMP=$OUTPUT"_temp"

if [ ! -d $TEMP ]; then mkdir -p $TEMP ; fi

# echo $SAMTOOLS
if [ ! -f $BAM.bai ];then 
  # echo "index bam...."
  $SAMTOOLS index $BAM 
fi

# echo "annotating..." 	
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
	rm $BAM
	echo "$msg End $0: `date`"
	exit 0
else
	echo "$msg End $0: `date`"
	exit 1
fi

echo "GATKannotationDONE" > $VCF.gatk.vcf.list

###-----check dir setting-----
if [ ! -d $cwd/log ]; then 
  echo -e $ERR"log dir not there!" 
  exit
fi
logDir=$cwd/log

##-----invoke merge VCFs
cnt=`ls -1 $cwd/split*gatk.vcf 2>/dev/null |wc -l`
if [[ -f $cwd/Moniter.txt ]] && [[ $cnt -gt 0 ]] ; then
  echo "Moniter file exisit!" 
else
  echo "Moniter var.vcf.gatk.vcf"  > $cwd/Moniter.txt
  echo "$srcDir/step-4_mergeVCF.sh $cwd $outputDir" | qsub -l mem=4g,time=140:: -e $logDir/ -o $logDir/ -cwd -N mergeVCF > $cwd/qsubMerge.log
  cat $tempDir/qsubMerge.log
fi 


echo "#--END---$(date)"
