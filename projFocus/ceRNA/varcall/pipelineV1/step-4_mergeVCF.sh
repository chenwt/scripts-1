#!/bin/bash
#$ -m e
#By: J.He
### input: working directory  
### output: <var.vcf.gatk.vcf> 
##Desp.: combine all subVCFs, must be executed in the working folder  
##input: temp-folder for the bam file; output directory full path

tempDir=$1
outputDir=$2
echo "Moniter var.vcf.gatk.vcf" > $tempDir/Moniter.txt
##---setting parameters-
sh /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh 
REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/humanGenome_hg19_regions104.bed
VCFCONCAT=~/tools/vcftools/vcftools_current/bin/vcf-concat
srcDIR=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
ERR="[ERR:]"
MSG="[MSG:]"
DONE="[DONE:]"
rootd=`pwd`
outVCFname=`echo $tempDir|awk -F"/" '{gsub("temp-","",$NF);print $NF} '`".var.vcf.gatk.vcf"

output=$outputDir/$outVCFname
echo $output
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
echo -e "total sub file:\t"$cnt
cntNF=`ls -1 split_*bam.var.vcf.gatk.vcf 2>/dev/null | wc -l`
cntjob=0
while [ $cntNF -le $cnt ] 
do
    if [ $cntNF -eq $cnt ] ; then
	 vcf_array=( `ls split_*bam.var.vcf.gatk.vcf ` )
	 echo ${vcf_array[@]}
	 ### merge VCFs 
	 $VCFCONCAT ${vcf_array[@]} > $output
	 break 
    else
      echo "$cntNF out of $cnt VCF files generated...$(date)"
    fi
    sleep 5m
    cntNF=`ls -1 split_*bam.var.vcf.gatk.vcf 2>/dev/null | wc -l`
done

if [[ $? == 0 ]] && [[ -f $output ]]
then
  echo -e "Deleting temp folder" 
  rm -R $tempDir
  echo "FINAL VCF GENERATED at $output"
  exit 0
else 
  echo -e "ERROR exit $(date)"
  exit 1
fi 
