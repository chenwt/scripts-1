#!/bin/bash
#$ -l mem=4g,time=24
#$ -N testRe
# input: bam file name 
# output:
# example: ~/scripts/projRe/do_02_cuff.sh TCGA-BH-A1ES-01A-11R-A137-07.bam
# echo "~/scripts/projRe/do_02_cuff.sh TCGA-BH-A1ES-01A-11R-A137-07.bam" | qsub -l mem=4g,time=8:: -cwd -N cufftest -e .log/ -e .log/

bam=$1
# cufflinks=$(which cufflinks)
cufflinks=/ifs/data/c2b2/ngs_lab/ngs/usr/bin/cufflinks
ind="/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/data/normal/"
gtf="/ifs/scratch/c2b2/ac_lab/jh3283/database/ncRNAlib/broad/lincRNAs_transcripts.gtf"
outrd="/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/"
outsd=$(echo ${bam} |sed -e "s/\.bam//g")
# echo ${cufflinks}
# echo ${gtf}
# echo $outsd

mkdir ${outrd}${outsd}
cmd="$cufflinks -G ${gtf} -o ${outrd}${outsd} ${ind}${bam}"
echo $cmd
$cmd