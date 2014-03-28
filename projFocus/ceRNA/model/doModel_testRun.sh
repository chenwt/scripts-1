#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 

#init
SUC="##--RUN-Suc"
ERR="##--RUN-Err"
CDT="Mar-7-2014"

##---testrun on different genes
# g='APC'
# g='BRAF'
# g='GATA6'
# g='NOTCH1'
g='CHEK1'
gDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/$g-temp
/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/getGeneData.sh $g $gDir


if [[ $? != 0 ]]
then
  echo $ERR
  exit
else
  echo $SUC 
fi
