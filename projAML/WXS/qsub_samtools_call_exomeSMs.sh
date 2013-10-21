#!/bin/bash
#Author: Jing He
#Date: Thu Jun 13 17:35:42 2013
#Last Update:  Thu Jun 13 17:35:01 UTC 2013
#Example:qsub_samtools_call_exomeSMs.sh <PID.txt> <outDir>
#        qsub_samtools_call_exomeSMs.sh <PID.txt> <outDir>
#NOTE: Should have a log directory under current working directory


 if [ $# >= 1 ]; then

     # scripts directory
     sd=$( dirname $( readlink -m $0 ) )

     res=$( readlink -m $2 )

     # loop through the cases--prepare a file listing all PID(with pairs)
     # for Tumor vs control, recurrent tumor vs control, should change the sample code below, tumor(04), recurrent tumor(09)
     while read PID
     do
       echo $PID
         if [ -e ${res}/jobs/callvar.${i}.job ]; then
             echo "${res}/jobs/callvar.${i}.job already exists"
         else
	   echo "${sd}/do_samtools_call_exomeSMs.sh $PID" | qsub -N callvar.${PID} -l mem=8G,time=20:: -e ${res}/log -o ${res}/log > ${res}/jobs/callvar.${PID}.job

           cat ${res}/jobs/callvar.${PID}.job
         fi

     done < $1

 else

     echo "error - need to supply an argument"

 fi
