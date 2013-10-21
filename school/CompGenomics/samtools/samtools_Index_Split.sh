#!/bin/bash
#Author: Jing He
#Date: Mar. 2013
#Last Update: Apr. 10, 2013
#Example: ~/scripts/school/compGenomics/samtools/run_samtools.sh
	###------------------------------indexing 
	# samtools index tumor.bam
	
	
    qsub -N indexT -e ./logs -o ./logs -b y -sync y -l mem=10G,time=10:: -S /bin/sh -cwd "samtools index tumor.bam"
    ###------------------------------splitting the bam files into different chromosomes
	
	for chrom in `seq 1 22` X Y
	do
		qsub -N splitT$chrom -e ./logs -o ./logs -b y -l mem=10G,time=12:: -S /bin/sh -cwd "samtools view -bh tumor.bam $chrom > tumor.$chrom.bam"
	done

