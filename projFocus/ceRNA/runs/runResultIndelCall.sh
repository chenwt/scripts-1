#!/bin/bash
#By: J.He
#TODO: 

###----test-----
#echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v2.sh /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/indelCall/test/test.bam"|qsub -l mem=8g,time=40:: -N indelTestrun -cwd -o ./log -e ./log >>qsub.log 
#tail -1 qsub.log

#-----test for one sample---------------------
##TCGA-A1-A0SD-01A-11D-A10Y-09.bam
#echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v2.sh $outDir /ifs/scratch/c2b2/TCGA/data/BRCA/WXS/TCGA-A1-A0SD-01A-11D-A10Y-09.bam "|qsub -l mem=8g,time=40:: -N indelTestrun -cwd -o ./log -e ./log >>qsub.log
#tail -1 qsub.log
#-------test--end----------------------------

##--making bamlist files--------------------------------------------
#cntTotalBam=`ls ${dataDir}/*bam|awk 'END{print NR}'`
#cnt=`echo "$cntTotalBam/20+1"|bc`
#cntPart=1
#echo -e "$cntTotalBam\t$cnt\t$cntPart"
#ls ${dataDir}/*bam > inputIndelCallBamList 
#for i in `seq 1 $cnt`
#do
# cntx=`echo $i*20|bc`
# cntm=`echo "$i*20-20"|bc`
# echo -e "line start:$cntm\tline end:$cntx"
# awk -v cntx=$cntx -v cntm=$cntm 'NR<=cntx&&NR>cntm{print $0}' inputIndelCallBamList > input.bam.list.${i}
#done
###----------------make--input-list-file--end--------------------------

#-------function-to-submit---bam-file-list---------
qsubRun (){
  cnt=0
  while read line
  do
    let cnt=$cnt+1
    jobid=`echo $line|awk 'BEGIN{FS="/"}{print substr($NF,9,4)}'`"_$cnt"
    echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v3.sh $outDir $line"|qsub -l mem=8g,time=160:: -N $jobid -cwd -o ./log -e ./log >>qsub.log
    tail -1 qsub.log
  done <$1
}

delBam() {
  bamlist=$1
  cnt=0
  for line in `ls *VCF.filtered.VCF`
  do
    bam=`echo $line|awk 'BEGIN{FS="_"}{print $1}'`
    pathBam=`grep $bam $bamlist`
    if [[ ! -z $pathBam ]]; then
      echo "deleting $pathBam"
      rm $pathBam 
      rm $pathBam.bai
      let cnt=$cnt+1
    fi
  done 

  echo -e "$cnt bam file deleted"
  echo "#----DONE----"
}

###-------function----end--------

#----data folder:
dataDir=/ifs/scratch/c2b2/TCGA/data/BRCA/WXS
outDir=$(pwd)

###---------run-----using--input---file---list
#qsubRun input.bam.list.7 
#qsubRun input.bam.list.0 
#qsubRun input.bam.list.1 
#qsubRun input.bam.list.2 
#qsubRun input.bam.list.3 
# qsubRun input.bam.list.brca_wxsInwgsDownloadlist2AND3.bam.tumor_02052014.txt_part1_10
# qsubRun input.bam.list.brca_wxsInwgsDownloadlist2AND3.bam.tumor_02052014.txt_part2_10
##------g

#delBam input.bam.list.0
#delBam input.bam.list.1
#delBam input.bam.list.2
delBam input.bam.list.3
delBam input.bam.list.brca_wxsInwgsDownloadlist2AND3.bam.tumor_02052014.txt_part1_10 

# touch input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt
# for inputlist in input.bam.list.4 input.bam.list.5 input.bam.list.6
# do
#   grep -f ../../data/sampleInfo/brca_wgs_bam_summary_02042014.tsv_TumorSample $inputlist >> input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt
# done 

# qsubRun input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt
# echo "/ifs/scratch/c2b2/TCGA/data/BRCA/WXS/TCGA-A7-A26G-01A-21D-A167-09.bam" > input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt_rescue1
# qsubRun input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt_rescue1


###------whole genome sequence
# cd ../../data/wgs
#grep -f brca_wgs_bam_summary_02042014.tsv_TumorSampleOverlapRNAseq /ifs/scratch/c2b2/TCGA/data/BRCA/WXS/input_gtdownload_v2.txt_part2|cut -f1 >> ../../result/indelCall/input.bam.list.brca_wxsInwgsDownloadlist2AND3.bam.tumor_02052014.txt
