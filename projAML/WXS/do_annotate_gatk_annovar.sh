#!/bin/bash
#$ -cwd
# J.HE
# TODO:gatk annotation requires more time and annovar requires more memory. So, better to seperate those two jobs
##Annotate the varaints files using GATK
VCF=$1
BAM1=$2
BAM2=$3
GLOBAL='/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/global_setting_GRC37lite.sh'
sc=$(pwd)
logd=${sc}/log
# VCF=`echo ${line} |awk '{print $1}'`
# BAMs=`echo ${line}|awk '{print $2,$3} '`

. ${GLOBAL}

echo $BPATH
CMD="sh $BPATH/gatk_annotate_variants_Jing.sh $VCF $BAM1 $BAM2 $GLOBAL "
echo $CMD
$BPATH/gatk_annotate_variants_Jing.sh $VCF $BAM1 $BAM2 $GLOBAL 
# $CMD
if [[ $? == 0 ]]
then
        echo "Success GATK-ANNOTATE"
        #Annotate the variants further using ANNOVAR genomic annotation
        #CMD="sh $BPATH/do_annovar_all.sh $VCF.gatk.vcf $VCF.gatk.vcf.annovar.vcf"            
        #echo $CMD
	jobname=`echo $VCF|awk 'BEGIN{FS="/"} {print $(NF-1)}'`
        echo "/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/do_annovar_all.sh $VCF.gatk.vcf $VCF.gatk.vcf.annovar.vcf" |qsub -l mem=20g,time=6:: -N annovar_${jobname} -e ${logd} -o ${logd} -cwd >> ${logd}/qsub.log
        tail -1 ${logd}/qsub.log
        # $CMD
        # $CMD 1> $VCF.log.annovar.o  2> $VCF.log.annovar.e
        
else
        echo "Failed GATK-ANNOTATE"
        echo "End $0: `date`"
        exit 1
fi
