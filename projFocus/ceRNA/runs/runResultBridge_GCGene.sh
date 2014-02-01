#!/bin/bash
#By: J.He
#TODO: 
#$ -cwd

##----------function----
function getGeneCerna(){
  ##----
  gene=$1
  outfile=$2
  head -1 $cernaNet > ${outfile}
  awk -v g=$gene 'BEGIN{FS=OFS="\t"}
		  NR>1&&($1==g||$2==g){
		print $0}' $cernaNet >> ${outfile}
}

function getCerna(){

  genelist=$1
  outfile=$2
  while read gene 
  do
    getGeneCerna $gene $outfile
  done < ${outfile}
  echo "$outfile ceRNA network extracted...."
}
##------function--end
wd=`pwd`
cernaNet=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/interactom/brca_ceRNA_network.txt
genelist=${wd}/gene.list
getGeneCerna CCND1 CCND1_ceRNANet.txt

