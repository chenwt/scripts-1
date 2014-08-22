#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

dataDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/data
dataOut=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfCerna/result/runAug21
cd $dataOut
if [ ! -d $dataOut/log ] ; then mkdir $dataOut/log ; fi
cnt=0
for f in `ls $dataDir/input*`
do
  tgene=`echo $f|awk -F"_" '{print $NF}'`
  echo $tgene
  # echo "~/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/tf/model/step8-5_reg_tfCerna.r --ctar $tgene --nboot 1000" 
  echo "~/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/tf/model/step8-5_reg_tfCerna.r --ctar $tgene --nboot 1000" | qsub -l mem=8g,time=48:: -e ./log -o ./log -N btReg_${tgene}  -cwd
  ((cnt++))
done

echo "job submmitte "$cnt
