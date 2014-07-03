#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:


##---- FOXC1 ceRNA network
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jul2014
cenet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_ceRNA_network.txt
lasnet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/candiReg/runApr30/summary/keyRegSummary_allRunMay_05212014_0.01
grednet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul1/optCorr.result_flex_max_1000.tsv.significant.summary.netfile

function extract {
  grep -w FOXC1 $cenet > $CWD/net.cernet
  
  grep -w "^FOXC1" $lasnet |awk '{split($2,a,";");for (i in a) {print $1"\t"a[i]}}'  >$CWD/net.lassonet
  
  grep -w "^FOXC1" $grednet |awk '{split($2,a,";");for (i in a) {print $1"\t"a[i]}}' >$CWD/net.greedynet
  
mutg=/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jul2014/FOX1.driverMut.gene
mutmibs=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/bs_with_addedBS/brca_som_startMut.miBS.sort.uniq.bed
muttfbs=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-1_bsSite/mut_in_TFBS/brca_som_start.tfBS.sorted.uni.mut.bed
mutall=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.bed

grep -wf $mutg $mutmibs > $CWD/FOX1.driverMut.gene.miBS.bed  
grep -wf $mutg $muttfbs > $CWD/FOX1.driverMut.gene.tfBS.bed  
grep -wf $mutg $mutall > $CWD/FOX1.driverMut.gene.all.bed  

}

# awk 'NR==FNR{a[$2]=$1;next;}{print a[$4]}' $CWD/net.lassonet $mutall |sort|uniq > $CWD/net.lassonet.mutated
awk 'FNR==NR{a[$2]=$2;next}{print a[$4]}' net.lassonet /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byPoint.matrix.bed |sort|uniq|sed 1d > net.lassonet.mutated



