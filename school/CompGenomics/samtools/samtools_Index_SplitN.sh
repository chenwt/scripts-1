#!/bin/bash
#Author: Jing He
#Date: Mar. 2013
#Last Update: Apr. 10, 2013
#Example: ~/scripts/school/compGenomics/samtools/run_samtools.sh
	###------------------------------indexing 
     
    qsub -N indexN -e ./logs -o ./logs -b y -sync y -l mem=10G,time=6:: -S /bin/sh -cwd "samtools index normal.bam"

	###------------------------------splitting the bam files into different chromosomes
	for chrom in `seq 1 22` X Y
	do
		qsub -N splitN$chrom -e ./logs -o ./logs -b y -l mem=8G,time=10:: -S /bin/sh -cwd  "samtools view -bh normal.bam $chrom > normal.$chrom.bam"
	done

