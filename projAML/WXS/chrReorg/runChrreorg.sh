#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 
pid="PARGVC"
/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/chrReorg/breakdancerPipeline.sh 
# echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/chrReorg/breakdancerPipeline.sh " |qsub -l mem=8g,time=20:: \
# -N brkD-${pid} -cwd -e temp-${pid} > logs.qsub
# tail -1 logs.qsub


