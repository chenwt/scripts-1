#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

RTRACE=/ifs/home/c2b2/ac_lab/rs3412/tools/TRUP/bin/RTrace.pl

samplefile=/ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/062014/tcgaRNAseq/fastq.wgs.ptumor.txt
READPOOL=/ifs/scratch/c2b2/TCGA/data/lincrna/rnaseq/brca/fastq/
AD=/ifs/home/c2b2/ac_lab/rs3412/tools/TRUP/annovar/
TH=2

for sample in `cat /ifs/data/c2b2/ac_lab/jh3283/projMisc/varReg/data/062014/tcgaRNAseq/fastq.wgs.ptumor.txt` 
  sample=`ls ${sample}*_1*|awk -F"/|_" '{print $(NF-1)}'`
  ###quality check
  /usr/bin/perl RTrace.pl --runlevel 1 --sampleName SAMPLE --readpool RP --root PD --threads TH --anno $AD 2>>run.log 

  ###mapping
  /usr/bin/perl RTrace.pl --runlevel 2 --sampleName SAMPLE --readpool RP --$root PD --anno $AD --threads TH --WIG --patient ID --tissue type --threads TH --gf pdf 2>>run.log


  ## gene/isoform quantification
  perl RTrace.pl --runlevel 5 --sampleName SAMPLE --readpool RP --root PD --anno $AD --threads TH --gtf-guide --known-trans refseq 2>>run.log

  perl $RTRACE --tmpDir /ifs/scratch/c2b2/ac_lab/rs3412/tmp/ --runlevel 1 --sampleName $SAMPLE --seqType p --root /ifs/scratch/shares/ngs-Califano/NET/RNA-seq/trup/100M_PE/ --anno $AD --readpool $READPOOL --gf pdf --threads 4 --Rbinary Rtrup --WIG 2>>/ifs/scratch/shares/ngs-Califano/NET/RNA-seq/trup/100M_PE/run.log.$SAMPLE
