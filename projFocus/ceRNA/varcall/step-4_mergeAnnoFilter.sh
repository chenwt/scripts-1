#!/bin/bash
#$ -m e
#By: J.He
### input: working directory  
### output: <var.vcf.gatk.vcf> 
##Desp.: combine all subVCFs, must be executed in the working folder  
##input: temp-folder for the bam file; output directory full path

tempDir=$1
outputDir=$2

##---setting parameters-
. /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh 
REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/humanGenome_hg19_regions104.bed
VCFCONCAT=~/tools/vcftools/vcftools_current/bin/vcf-concat
srcDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
finalDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered
PYTHON=~/tools/python/Python_current/python
err="[ERR:]"
rootd=`pwd`
outVCFname=`echo $tempDir|awk -F"/" '{gsub("temp-","",$NF);print $NF} '`".var.vcf.gatk.vcf"
output=$outputDir/$outVCFname

###-----check dir setting-----
if [ ! -d $tempDir ]; then 
  echo -e $ERR"temp dir not there!" 
  exit
fi
  
if [ ! -d $tempDir/log ]; then 
  echo -e $ERR"log dir not there!" 
  exit
fi

cwd=$tempDir
logDir=$cwd/log


###--------check existing var.vcf.gatk.vcf----
cd $cwd
cnt=`awk 'END{print NR}' $REGIONFILE`
((cnt -=1 ))
echo -e "total sub file:\t"$cnt
cntNF=`ls -1 split_*bam.var.vcf.gatk.vcf 2>/dev/null | wc -l`
cntjob=0
#--
if [[ $cntNF -eq $cnt ]]; then
    $PYTHON $srcDIR/vcfConcat_local.py -i var.vcf.gatk.vcf$ -o $output
else
   echo "$err: vcf files not ready"
   exit
fi

if [[ $? == 0 ]] && [[ -f $output ]]
then
   cd $rootd
   echo -e "Deleting temp folder" 
   rm -R $tempDir 
   echo "FINAL VCF GENERATED at $output"
else 
  echo -e "$err vcf concating $(date)"
  exit 1
fi 

##--annovar
$srcDIR/step-5_doAnnovarAll.sh $output  

if [[ $? !=0 ]] && [[ ! -f $output.annovar.summary.genome_summary.csv.vcf ]] ; then
  echo "$err: annovar annotating error"
  exit 1 
fi 
rm $output 
rm $output.annovar.summary.exome_summary.csv
rm $output.annovar.summary.genome_summary.csv



##---fillter
outputFiltered=$output.summary.genome.csv.vcf.localfiltered.vcf
~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/filterMutation_local.py -i $output.summary.genome.csv.vcf -o $outputFiltered

filename=`echo $outputFiltered |awk -F"/" '{print $NF}'`

if [[ $? !=0 ]] && [[ ! -f $fileName ]] ; then
  echo "ERR: final step cleaning files" 
fi

##---copy and clean file

cp $fileName $finalDir/

if [[ $? !=0 ]] && [[ ! -f $fileName ]] ; then
  echo "ERR: final step cleaning files" 
fi
rm $fileName 
echo -e "##-------FINAL OUTPUT $finalDir/$fileName "
