#!/bin/bash
# J.HE
#$ -cwd
#$ -l mem=4g,time=10::

gatk="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar"
ref="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa"
# input=
bamd="/ifs/scratch/c2b2/TCGA/data/BRCA/WGS/" 
vcfd=
pid=$1
chr=$1
input1=$(ls ${bamd}*${pid}*_split/${chr}.bam | awk '{print $1}') 
input2=$(ls ${bamd}*${pid}*_split/${chr}.bam | awk '{print $2}') 
# input3=
wd=
cd wd
cmd="java -Xmx2g -jar $gatk -R $ref \\
 -V ${output_folder}/$s.$chr.vcf  \\
 -T VariantAnnotator  \\
 -L ${chr} \\ 
 -I ${input1} \\ 
 -I ${input2} \\
 -I ${input3} \\
 -A DepthPerAlleleBySample 
 -o ${output_folder}/$s.$chr.gatk.vcf"
 echo $cmd