#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:




##---get input ceRNA target genes
tfmutDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/test
keyRegfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summary/tgeneKeyreg_0.01

awk 'NR>1&&NR<=11{print $1}'  ${keyRegfile} > ${tfmutDir}/input.tfmut.test 




