#!/bin/bash
#$ -cwd
#$ -l mem=8G,time=20::
#$ -o ./log -j yes

#Author:J.He
#while read pid 
#do
# if running by recurrent tumor versus normal
#./samtools_call_exomeSMs.sh ${pid}-09A ${pid}-14A

# if running by tumor versus normal
#./samtools_call_exomeSMs.sh ${pid}-04A ${pid}-14A

# if running by recurrent tumor versus tumor
#./samtools_call_exomeSMs.sh ${pid}-09A ${pid}-04A
#done < PID



./samtools_call_exomeSMs.sh PAEEYP-09A PAEEYP-04A
