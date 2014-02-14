#!/bin/bash
#By: J.He
#$ -cwd

##-D output per-sample read depth -S per-sample Phre-scaled strand bias p-value -u for piping -f faidx index reference file 
sh /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh
SAMTOOLS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools
BCFTOOLS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bcftools
VCFUTILS=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/vcfutils.pl
$SAMTOOLS mpileup -DSuf $REF $1  |$BCFTOOLS  view -bvcg -p 0.9 - > $1.var.bcf
$BCFTOOLS view $1.var.bcf | $VCFUTILS varFilter -D 200 > $1.var.vcf
rm $1.var.bcf
echo $1.var.vcf > $1.var.vcf.list
