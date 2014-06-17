#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

PYTHON=~/tools/python/Python_current/python
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som

firstStep() {
  tcgamaf=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/03102014/tcgal2som/genome.wustl.edu__Illumina_All.maf
  # awk -F"\t" 'NR>1{print $16}' $tcgamaf |sort|uniq  > $CWD/all_Tumor_barcode.list
  out=$CWD/brca_somlevel2_byPoint.matrix
  # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/genSomMatFromTCGAl2_v2.py -i $tcgamaf -o $out 
  out=$CWD/brca_somlevel2_byGene.matrix
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/step1-1_getSomMutMatrix.py -i $tcgamaf -o $out -t gene 
}

##------extract only smps which have expression data

# head -1 brca_somlevel2_byGene.matrix |tr "\t" "\n" > barcode_allMut.list
# grep -f ../gslist/barcode_tumor_Mar-21-2014.list barcode_allMut.list > barcode_allMut_in_exp.list
# sed -i "1igene" barcode_allMut_in_exp.list
# fullMut=$CWD/brca_somlevel2_byGene.matrix
# ~/bin/trfile $fullMut $fullMut.temp
# ~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/extractTargetMatFromAllMat.py -i barcode_allMut_in_exp.list -d $fullMut.temp -o $fullMut.inExpSmp.matrix.temp
# ~/bin/trfile $fullMut.inExpSmp.matrix.temp $fullMut.inExpSmp.matrix
# rm *temp
# head -1 brca_somlevel2_byGene.matrix.inExpSmp.matrix > header_inExpSmp
# cat brca_somlevel2_byGene.matrix.inExpSmp.matrix | awk '{sum=0;for(i=2;i<=NF;i++)sum+=$i;if(sum>0||NR==1) print $0}' > brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero

##---------summary
# cat brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero |awk '{sum=0;for(i=2;i<=NF;i++)sum+=$i;if(sum>0||NR==1) print $1"\t"sum}' > brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero.geneMutfreq
# ~/bin/trfile brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero.trfile
# cat brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero.trfile | awk '{sum=0;for(i=2;i<=NF;i++)sum+=$i;if(sum>0||NR==1) print $1"\t"sum}' > brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero.smpMutfreq
# rm *trfile *temp
# sort -k 2nr  brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero.geneMutfreq | cut -f2 |sort|uniq -c |sort -k 1nr|less

###-----get promoter region mutation, 3putr mutation
function getRegionSpecMut {
  CDT=`date|awk '{print $2"-"$3"-"$6}'`
  input=$CWD/brca_somlevel2_byPoint.matrix
  tssfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar20_Tss.tsv.signle.tsv
  utr3pfile=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar20_3pUTR.tsv

  promotersize=2000
  output=$input.promoter2k.$CDT.matrix
  cmd="$PYTHON  $crnsc/processData/selectTargetRegionRow.py -i $input -t $tssfile -c $promotersize -o $output &"
  echo $cmd
  $cmd 

  promotersize=1000
  output=$input.promoter1k.$CDT.matrix
  cmd="$PYTHON $crnsc/processData/selectTargetRegionRow.py -i $input -t $tssfile -c $promotersize -o $output &" 
  echo $cmd
  $cmd

  size=0
  output=$input.utr3p.$CDT.matrix
  cmd="$PYTHON $crnsc/processData/selectTargetRegionRow.py -i $input -t $utr3pfile -c $size -o $output &"
  echo $cmd
  $cmd

}

getRegionSpecMut
