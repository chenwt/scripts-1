#!/bin/bash
#$ -l mem=4g,time=48::
#$ -o ./log -e ./log
#$ -cwd
#Usage:
#Example
#Author:J.He
#Date: Jun 20, 2013


while read pid
do 
#sd=$(dirname $0)
sd=$(pwd)

###--------------------------------------
## for 04-09-14 style
cmd="sh ${sd}/runMuTect.sh ${pid}-04A ${pid}-14A"
echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-04a-14a -cwd 

cmd="sh ${sd}/runMuTect.sh ${pid}-09A ${pid}-14A"
echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-09a-14a -cwd 
#
#cmd="sh ${sd}/runMuTect.sh ${pid}-09A ${pid}-04A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-09a-04a -cwd 

#rerun this section
cmd="sh ${sd}/runMuTect.sh ${pid}-04A ${pid}-09A"
echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-04a-09a -cwd 




###--------------------------------------
## for 03 04 14 
#cmd="sh ${sd}/runMuTect.sh ${pid}-04A ${pid}-14A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-04a-14a -cwd 
#
#cmd="sh ${sd}/runMuTect.sh ${pid}-03A ${pid}-14A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-09a-14a -cwd 
#
#cmd="sh ${sd}/runMuTect.sh ${pid}-04A ${pid}-03A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-09a-04a -cwd 
#

###--------------------------------------
## for 09 40 14 
#cmd="sh ${sd}/runMuTect.sh ${pid}-09A ${pid}-14A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-09a-14a -cwd 
##
#cmd="sh ${sd}/runMuTect.sh ${pid}-40A ${pid}-14A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-40a-14a -cwd 
##
#cmd="sh ${sd}/runMuTect.sh ${pid}-40A ${pid}-09A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-40a-09a -cwd 
#


###--------------------------------------
## for 03 04 10 
#cmd="sh ${sd}/runMuTect.sh ${pid}-03A ${pid}-10A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-03a-10a -cwd 
##
#cmd="sh ${sd}/runMuTect.sh ${pid}-04A ${pid}-10A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-04a-10a -cwd 
##
#cmd="sh ${sd}/runMuTect.sh ${pid}-04A ${pid}-03A"
#echo $cmd| qsub -l mem=4g,time=48:: -o ./log -e ./log -N mutect_${pid}-04a-03a -cwd 
##


done < PID_input.txt 
