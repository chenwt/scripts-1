#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:
source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step4_MSMC/weighted

keyRegDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/sigKeyReg
gslistf=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list
mutfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix
outputf=$CWD/brca_mutCernaDriver_wgsc_${CDT}
# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_WGSC.py -a 0.8 -t mean -d $keyRegDir -g $gslistf -m $mutfile -o $outputf &
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_WGSC.py -a 0.8 -t max -d $keyRegDir -g $gslistf -m $mutfile -o $outputf &
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_WGSC.py -a 0.85 -t mean -d $keyRegDir -g $gslistf -m $mutfile -o $outputf &
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step3-4_WGSC.py -a 0.85 -t max -d $keyRegDir -g $gslistf -m $mutfile -o $outputf &

