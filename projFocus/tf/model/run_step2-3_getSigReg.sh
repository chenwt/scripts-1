#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

PYTHON=~/tools/python/Python_current/python
mysc=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/tf/model/step2-3_extract_sigReg_from_candiReg.py

for file in `ls $1/*_candidateRegs.txt`
do 
  $PYTHON $mysc $file
done 

# echo "[END]"
