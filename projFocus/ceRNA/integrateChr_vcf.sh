#!/bin/bash
#$ -cwd
#$ -l mem=4g,time=2::
# J.HE
# Integrate samtools chr wised calling result
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/
sc=$(basename $0)
echo "Start $0: `date` " >> ${wd}log/log_${sc}
pid=$1
# while read pid
# do 
 	# cmd="grep ^# ${wd}${pid}/out_pairY.vcf > ${wd}${pid}/header.vcf"
 	grep ^\# ${wd}${pid}/out_pairY.vcf > ${wd}${pid}/header.vcf
 	echo $cmd >> ${wd}log/log_${sc}
 	$cmd
 	if [ -f ${wd}${pid}/temp.joint.vcf ]; then 
 		rm ${wd}${pid}/temp.joint.vcf 
 	fi
	for chr in $(seq 1 22) X Y
	do
		cmd="grep -v ^# ${wd}${pid}/out_pair${chr}.vcf > ${wd}${pid}/temp.joint.vcf"
		echo $cmd >> ${wd}/log/log_${sc}
		# $cmd
 		grep -v ^\# ${wd}${pid}/out_pair${chr}.vcf >> ${wd}${pid}/temp.joint.vcf
	done
	cmd="cat ${wd}${pid}/header.vcf ${wd}${pid}/temp.joint.vcf > ${wd}${pid}.vcf"
	echo $cmd >> ${wd}log/log_${sc}
 	$cmd
 	echo "$pid.vcf generated!"
# done < ${wd}input_integrateChr_vcf.txt
echo "End `date` " >> ${wd}log/log_${sc}