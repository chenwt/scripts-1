#!/bin/bash
#$ -l mem=4g,time=24
#$ -N testRe
# input ： no
# require: <file os bam filename list,one per line M rows>, <bed file N rows> , files on server not local
# output:two files: <tab seperated txt file(N * (M+1) matrix)>
# USAGE: sh get_tssCount.sh 
# TODO: ADD control to generate big file, need moniter the temp_awk* file number

mybed=/ifs/scratch/c2b2/ac_lab/jh3283/database/ncRNAlib/broad/lincRNAs_transcripts.bed
wd=/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/data/
tempd=/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/temp/
logd=/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/log/
resd=/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/
cnt=1
total=$(ls *bam|wc -l)
# for bam in $(ls *bam)
# do 
# 	cmd="~/tools/bedtools/multiBamCov -p -bams ${bam} -bed ${mybed}  >  ${tempd}temp_${cnt}"
# 	if [ $cnt == $total ]; then
# 	   	# echo "last job"
# 	   	echo $cmd | qsub -l mem=4g,time=2:: -e ${logd} -o $logd -N count${cnt} -sync y -cwd >> ${logd}log_qsub
# 		tail -1 ${logd}log_qsub 
# 	else 
# 	   	echo $cmd | qsub -l mem=4g,time=2:: -e ${logd} -o $logd -N count${cnt} -cwd >> ${logd}log_qsub
# 	   	tail -1 ${logd}log_qsub 
# 	   	# echo "not the last"
# 	fi
#    	cnt=$((cnt + 1))
# done

# for i in $(seq 1 ${total})
# do
# 	if [ ${i} == 1 ]; then
# 		awk 'BEGIN{OFS="\t"}{print $1"_"$2"_"$3,$NF}' ${tempd}temp_${i} > ${tempd}temp_awk${i}
# 	else 
# 		awk 'BEGIN{OFS="\t"}{print $NF}' ${tempd}temp_${i}> ${tempd}temp_awk${i}
# 	fi
# done


# echo "Done! need you to add header!!"

paste ${tempd}temp_awk* > ${resd}result_tss_count_${total}.txt

### not need any more
# ls *bam| tr "*\n" "\t" | awk 'BEGin{OFS="\t"}{print "CHR","START","END",$0}'
# echo "chr	start	end	$fns"
# echo "CHR, START, END, $fns" testfile.csv