#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist
gslistDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist

source /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/geneUtilsRuns.sh 
tfnet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_geneid-regulon.rda.4col.txt.naomit.symbol
tfnet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist/brca_tcga_rnaseq851_aracne.net.gata3

# tglist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist/brca_tcga_rnaseq851_aracne.tftarget.list

tglist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist/brca_tcga_rnaseq851_aracne.tftarget.list.gata3

function makeGslist {
  awk 'NR>1{print $1}' $tfnet | sort|uniq > $gslistDir/brca_tcga_rnaseq851_aracne.tflist
  expfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
  
  awk 'NR==1' $expfile |tr "\t" "; " > $gslistDir/brca_tcga_allsample.idlist
  ~/tools/python/Python_current/python makeGslist.py $gslist
}


###----change to tf targets
function makeGslistTarget {
  awk 'NR>1{print $2}' $tfnet | sort|uniq > $tglist 
  # expfile=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix
  # awk 'NR==2' $expfile |tr "\t" "; " > $CWD/brca_tcga_allsample.idlist
  ~/tools/python/Python_current/python makeGslist.py $tglist
}

# makeGslistTarget 

### after runnning all step1, get all matrix ready, do this script
# tglist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist/brca_tcga_rnaseq851_aracne.tflist


###-------------Dysregulated (no need)

# # gslist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/gslist/brca_tcga_rnaseq851_aracne.tflistallsample.gslist

gslist=${tglist}_allsample.gslist
echo $gslist

function missGenomic {
  tumor=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix
  normal=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix
  out=${gslist}.deg
  $RSCRIPT /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-2_getDEgeneSample.r --tumor $tumor --normal $normal --gslist $gslist --out $out 
}

gslist=
awk 'NR==FNR{a[$1]=$0;next}{print a[$1]}' ${gslist} /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix|sort |uniq |sed 1d > ${gslist}.hasExp
gslist=${gslist}.hasExp 

# $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genStat4Gslist.py -i $gslist -c $tfnet -o ${gslist}_stat


tfnet=brca_tcga_rnaseq851_aracne.tftarget.list_allsample.gslist
###
$PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpCernet.py $tfnet

echo $gslist
echo "[END]"
