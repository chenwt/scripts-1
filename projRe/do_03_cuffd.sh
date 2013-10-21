#!/bin/bash
#$ -l mem=4g,time=4::
#$ -e ../log/ -o ../log/
#$ -N cuffd -cwd
# input: <gtf file from do_02_cuff.sh > <bam file name> 
# output:

# cuffmerge
gtf="/ifs/scratch/c2b2/ac_lab/jh3283/database/ncRNAlib/broad/lincRNAs_transcripts.gtf"
bamN='TCGA-A7-A0DB-11A-33R-A089-07.bam,TCGA-AC-A2FF-11A-13R-A17B-07.bam,TCGA-BH-A0BJ-11A-23R-A089-07.bam,TCGA-BH-A0BQ-11A-33R-A115-07.bam,TCGA-BH-A0DH-11A-31R-A089-07.bam'
bamT='TCGA-BH-A18V-01A-11R-A12D-07.bam,TCGA-BH-A1ES-01A-11R-A137-07.bam,TCGA-BH-A1FE-01A-11R-A13Q-07.bam,TCGA-E2-A15A-01A-11R-A12D-07.bam,TCGA-E2-A15E-01A-11R-A12D-07.bam,TCGA-E2-A15K-01A-11R-A12P-07.bam'
outd="/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/cuffd/"


cmd="cuffdiff -p 4 -o ${outd} ${gtf} ${bamT} ${bamC} "
echo ${cmd}
${cmd}