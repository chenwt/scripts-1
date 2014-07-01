#!/bin/bash
#$ -cwd
#By: J.He
#Desp.: mirBS data folder processing 
#Purpose: run
#last modified: 

source  $crnsc/geneUtilsRuns.sh
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS 

src=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/csv2bed_dataMirBS_mc2014.py

function localtest {
  head -1000 $CWD/MC2014_human_mirBS.bed.csv > $CWD/test.input
  $PYTHON $src $CWD/test.input $CWD/test.out.bed
  rm *test*
}

function getBed {
  input=$CWD/MC2014_human_mirBS.bed.csv 
  output=$CWD/miBS_mc2014_human.bed
  $PYTHON $src $input $output
}
# getBed


src=~/scripts/projFocus/ceRNA/processData/parclip2bed_dataMirBS_17319CCRs.py
function localtest {
 head -1000 $CWD/PAR-CLIP_17319CCRs  > $CWD/test.input
 $PYTHON $src $CWD/test.input $CWD/test.output
}
function getBed2 {
 $PYTHON $src $CWD/PAR-CLIP_17319CCRs $CWD/miBS_PARCLIP_17319ccrs.bed
}
# getBed2


###---same strategy as getting cupid binding site, do to get 3p utr binding site for CLASH data 
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/utrBS_to_genmoicsCoord/clash
function getClashBED {
  src=~/scripts/projFocus/ceRNA/processData/clash2bed_dataMirBS_6102sites.py
  rnaseq=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/3PrimeUTR
  # dnaseq=$CWD/refseq_human_hg19_DNA_3UTRexon_CLASH_geneName.fasta
  dnaseq=$CWD/refseq_human_hg19_clashgene_wholegene.fasta
  # $PYTHON $src $CWD/CLASH_6102sites $rnaseq $dnaseq $CWD/miBS_CLASH_6012sites_in_3primeUTR.bed
  $PYTHON $src $CWD/CLASH_6102sites $rnaseq $dnaseq $CWD/miBS_CLASH_6012sites.bed
  sort -k 1 -k 2n -k 3n $CWD/miBS_CLASH_6012sites.bed |uniq  > $CWD/miBS_CLASH_6012sites.uniq.bed
  cp $CWD/miBS_CLASH_6012sites.uniq.bed ../miBS_CLASH_6012sites.uniq.bed
  rm $CWD/miBS_CLASH_6012sites.bed
}

# getClashBED &


###----cupid again

# mkdir cupid
CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS/utrBS_to_genmoicsCoord/cupid

function splitData {
  # ~/bin/splitByN CupidPred.scored.geneName 3500
  
  # awk '{print $2"\t"$0}' $CWD/CupidPred.scored > $CWD/temp
  for i in `seq 1 5`
  do 
    ~/bin/grepf2f $CWD/CupidPred.scored.geneName_${i} $CWD/temp $CWD/CupidPred.scored.temp_${i} & 
  done
  for i in `seq 1 5`
  do
    awk '{print $2,$3,$4,$5,$6}' $CWD/CupidPred.scored.temp_${i}  > $CWD/CupidPred.scored_${i}
  done 
}

src=~/scripts/projFocus/ceRNA/processData/cupid2bed_dataMirBS_cupidScored.py

function getBed {
#potentially take longer time
    # i=1
    # input=$CWD/CupidPred.scored_${i}
    # rnaseq=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/3PrimeUTR
    # # dnaseq=$CWD/refseq_human_hg19_cupidGenename30bpupUTR3_${i}
    # dnaseq=$CWD/refseq_human_hg19_cupide_wholegene_${i}.fasta
    # $PYTHON $src $input $rnaseq $dnaseq ${input}.bed & 
  for i in `seq 3 5`
  do
    input=$CWD/CupidPred.scored_${i}
    rnaseq=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cupidPrediction/3PrimeUTR
    # dnaseq=$CWD/refseq_human_hg19_cupidGenename30bpupUTR3_${i}
    dnaseq=$CWD/refseq_human_hg19_cupide_wholegene_${i}.fasta
    $PYTHON $src $input $rnaseq $dnaseq ${input}.bed  &
  done 

}
# getBed 

function getUniqBS {
  for i in `seq 1 5`
  do
     sort -k 1 -k 2n -k 3n $CWD/CupidPred.scored_${i}.bed | cut -f 1-6| uniq |awk 'BEGIN{FS=OFS="\t"}
       {
        if (a[$1"\t"$2"\t"$3"\t"$4"\t"$5] < $6) a[$1"\t"$2"\t"$3"\t"$4"\t"$5]=$6;
       }
        END{for(i in a) {print i,a[i];}}' > $CWD/CupidPred.scored_${i}.bed.uniq 
  done 
}

# getUniqBS

function getUniqBS {
  for i in `seq 1 5`
  do
    cat $CWD/CupidPred.scored_${i}.bed.uniq >> $CWD/CupidPred.scored_all.bed
  done 
  
  sort -k 1 -k 2n -k 3n $CWD/CupidPred.scored_all.bed | uniq |awk 'BEGIN{FS=OFS="\t"}
         {
          if (a[$1"\t"$2"\t"$3"\t"$4"\t"$5] < $6) a[$1"\t"$2"\t"$3"\t"$4"\t"$5]=$6;
         }
          END{for(i in a) {print i,a[i];}}'|sort -k 1 -k 2n -k 3n |uniq  > $CWD/CupidPred.scored_all.bed.uniq 
  cp $CWD/CupidPred.scored_all.bed.uniq ../miBS_CupidScore.uniq.bed
}

# getUniqBS

CWD=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/miBS

function combAllBS {
  src=~/scripts/projFocus/ceRNA/processData/combineBS_bed.py
  clash=$CWD/miBS_CLASH_6012sites.uniq.bed
  parclip=$CWD/miBS_PARCLIP_17319ccrs.bed
  mc2004=$CWD/miBS_mc2014_human.bed
  # cupid=$CWD/miBS_CupidScore.uniq.bed
  cupid=$CWD/miBS_cupid.all.scored.bed
  output=$CWD/miBS_clash_parclip_mc2014_cupid.bed
  $PYTHON  $src $clash $parclip $mc2004 $cupid $output & 
  sort -k 1 -k 2n -k 3n $CWD/miBS_clash_parclip_mc2014_cupid.bed |uniq > $CWD/miBS_clash_parclip_mc2014_cupid.uniq.bed
}

combAllBS
