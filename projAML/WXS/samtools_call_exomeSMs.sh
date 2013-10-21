#!/bin/bash
#$ -cwd

#Author:J.He
#Usage: samtools_call_exomSMs.sh  <PID-SampleCode Tumor> <PID-SampleCode Control>
#Examples: samtools_call_exomSMs.sh PAYYEP-04 PAYYEP-14
#Description: this is called by do_samtools_calvar.sh, make sure all bam files are in the same working folder as the scripts 
##--------------------------------------
#environment variables
MYREF='/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa'
EXOBED='/ifs/scratch/c2b2/ac_lab/jh3283/ref/Exome_Targeted_Regions.BED'
BAM1=$1.rmdup.new.bam
BAM2=$2.rmdup.new.bam
echo $MYREF
echo $EXOBED
echo $BAM1"    "$BAM2


##--------------------------------------
#check BAM file
if [[ ! -f $1.rmdup.new.bam || ! -f $2.rmdup.new.bam ]] ;then  
  echo "NO BAM files provided!" 
  exit 1 
fi

###------------------------------call raw variants
samtools mpileup -ugSDf $MYREF -l $EXOBED $BAM1 $BAM2  | bcftools view -bvcgT pair -p 1.1 - > $1.var.bcf
echo $date  >> log$1
echo "pileup done"  >> log$1

bcftools view $1.var.bcf | vcfutils.pl varFilter -D 10000 > $1.var.vcf
echo "Calling $1 done!" >> log$1

##call somatic variantsq
ruby ~/scripts/school/compGenomics/samtools/do_filter-somatic.rb -v $1.var.vcf
echo "Filtering somatic variants $1 DONE!" >> log$1


#use annovar
sh $ANNOVAR/do_annovar_all.sh $1.var.vcf.somaticfiltered.vcf $1.var.vcf.somaticfiltered.annotated.vcf

echo "Annotation of $1 FINISHed!" > log$1
