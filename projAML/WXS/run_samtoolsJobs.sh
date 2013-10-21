#!/bin/bash
#Author: Jing He
#Date: Apr 16, 2013
#Last Update: Jun 10, 2013
#Example:  run_samtoolsJobs.sh PID

BASEDIR=$(dirname $(readlink -f $0))
 while read LINE
 do
	echo $(date) >> logsrun
	cd $BASEDIR/$LINE
	echo $LINE >> logsrun
#	 cp $BASEDIR/samtools_Index_Split.sh .
#	 cp $BASEDIR/samtools_Index_SplitN.sh . 
    qsub -N call$LINE$chrom -e ./logs -o ./logs -l mem=10G,time=24:: -S /bin/sh -cwd -v PID=$LINE samtools_Callvar_TN.sh
	done
#	
