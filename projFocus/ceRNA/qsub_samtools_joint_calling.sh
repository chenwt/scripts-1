#!/bin/bash
#$ -cwd
#$ -l mem=4g,time=48::
# input: <full path of working directory where the bam.list is>
# TODO: need to change the working directory 
inp=$1
GLOBAL='/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/global_setting_GRC37lite.sh'
# logd='/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/log'
if [ ! -d ${inp} ]; then
	echo "Directory for bam not exist!"
	exit 1
fi
cd ${inp}

if [ ! -d log ]; then
	mkdir log
fi
logd=$(readlink -m ./log)

#for i in $(seq 1 22) X Y 
for i in 1 2 
do 
	cmd="/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/samtools_joint_calling.sh \
		-I ${inp}/bams${i}.list \
		-g ${GLOBAL} \
		-O ${inp}/out_pair${i}.vcf \
		-m 8 \
		-A no"
	# echo ${cmd}
	# ${cmd}
	echo ${cmd} | qsub -l mem=8g,time=140:: -N call${i} -cwd -e ${logd}/ -o ${logd}/ >> ${logd}/log_qsub_call
	tail -1 ${logd}/log_qsub_call
done
