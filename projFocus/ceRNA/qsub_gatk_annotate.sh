#!/bin/bash
#

inp=$1
cd ${inp}
echo $(pwd)

if [ ! -d log ]; then
	echo " making log dir"
	mkdir log
fi
logd=$(readlink -m ./log)
echo ${logd}
for chr in $(seq 1 22) X Y
do 
	cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/gatk_annotate_variants.sh \
	${inp}/out_pair${chr}.vcf  \
	${inp}/bams${chr}.list \
	/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/global_setting_GRC37lite.sh"
	# echo ${cmd}
	echo ${cmd} | qsub -l mem=4g,time=10:: -N annot${chr} -e ${logd}/ -o ${logd}/  >> ${logd}/qsub_log
	tail -1  ${logd}/qsub_log
done 
