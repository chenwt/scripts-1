#!/bin/bash
#By: J.He
#TODO: 


##---testrun
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/test/filterIndel.py -i /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/TCGA-A1-A0SD-01A-11D-A10Y-09.bam_dindel_ouput.variantCalls.VCF -o /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/indelFinal/TCGA-A1-A0SD-01A-11D-A10Y-09.bam_dindel_ouput.variantCalls.VCF.filtered.VCF

PYTHON=~/tools/python/Python_current/python
inputDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall
outputDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/indelFinal


###----
for vcf in `ls $inputDir/*VCF`
do
  filename=`echo $vcf|awk 'BEGIN{FS="/"}{print $NF}'`
  echo "processing file..."
  $PYTHON ~/scripts/projFocus/ceRNA/filterIndel.py -i $vcf -o $outputDir/$filename.filtered.VCF
done 
echo -e "#---END---"
