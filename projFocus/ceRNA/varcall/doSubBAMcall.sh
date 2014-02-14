#!/bin/bash
#By: J.He
#$ -cwd 
#TODO: 
##input: <full path of bam > <bed file include the region to split the bam>
##output:

BAM=$1
##----split bam file
GLOBAL=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/global_setting_projFocueCeRNA.sh
. $GLOBAL
GATKJAR=/ifs/home/c2b2/ac_lab/jh3283/tools/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar
REGIONFILE=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/wgs/humanGenome_hg19_regions104.bed
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
bamName=`echo $BAM|awk 'BEGIN{FS="/"}{print $NF}'`
if [ ! -d temp-$bamName ]; then mkdir temp-$bamName ; fi
rootd=`pwd`
cwd=$rootd/temp-$bamName
if [ ! -d $cwd/log ]; then mkdir $cwd/log ; fi
logDir=$cwd/log
###--need 24 hours for run and 24 for waiting
/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall/subBAMCall.sh $REGIONFILE $BAM $cwd

echo "All subBAM calling job submitted!"
echo "#---DONE---"
