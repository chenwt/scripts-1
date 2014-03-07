#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
filtDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered

function annovarAll(){
#-----------check current dir, and redo annovar annotation for those which didn't do
  rawDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars
  cd $rawDir
  rawVCF_array=( `ls TCGA-*var.vcf.gatk.vcf` )
  echo ${rawVCF_array[@]}
  for tempVCF in "${rawVCF_array[@]}"
  do
      cntFiltVCF=`ls $tempVCF.annovar.summary.genome_summary.csv.vcf 2>/dev/null |wc -l`
      if [ $cntFiltVCF  == 0 ]
      then
      	  # echo "annovar annotating... $tempVCF"
	  tempVCFp=`readlink -f $tempVCF`
    	  pid=`echo $tempVCF|awk -F"-" '{print $3}' `
	  echo $pid
	  echo "$srcDir/step-5_doAnnovarAll.sh $tempVCFp " | qsub -l mem=14g,time=20:: -N af_$pid -e ./log -o ./log -cwd >>qsub.log
	  tail -1 qsub.log
      fi
  done
}

function annovarVCF(){
 pid=`echo $1|awk -F"-" '{print $3}' `
  echo "$srcDir/step-5_doAnnovarAll.sh $1 " | qsub -l mem=8g,time=10:: -N af_$pid -e ./log -o ./log -cwd >>qsub.log
  tail -1 qsub.log
}

function filterAll(){
  rawDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars
  filtDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered
  cd $rawDir
  rawVCF_array=( `ls TCGA-*var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf` )
  cnt=0
  for tempVCF in "${rawVCF_array[@]}"
  do
      temp=`echo $tempVCF|sed s/.annovar.summary.genome_summary.csv.vcf//g`
      # if [ -f $temp ]; then 
	# echo " deleting $temp" 
	# rm $temp 
      # fi
      cntFiltVCF=`ls $filtDir/$temp.annovar.summary.genome_summary.csv.vcf.filtered.vcf 2>/dev/null |wc -l`
      if [ $cntFiltVCF  == 0 ]
      then
    	  pid=`echo $tempVCF|awk -F"-" '{print $3}' `
	  (( cnt+=1 ))
	  echo " /ifs/home/c2b2/ac_lab/jh3283/tools/python/Python-2.7.5/python $srcDir/filterMutation_local.py -i $tempVCF -o $filtDir/$tempVCF.filtered.vcf" | qsub -l mem=10g,time=8:: -N filt_$pid -e $filtDir/log -o $filtDir/log -cwd >> $filtDir/qsub.log 
	  tail -1 $filtDir/qsub.log 
      fi
  done
  echo "filter job $cnt VCF"
}

checkVcf2Download(){
   vcfDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered
   downloadFile=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/wgs/input_gtBatch_v2.txt
   output=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/wgs/input_gtBatch_v2.txt_failedMar01.txt
   awk -F"\t" '{print $1}' $downloadFile > $downloadedFile.temp 
   ls $vcfDir/TCGA*vcf |awk -F"/" '{print $NF}'|sed s/.bam.var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf.filtered.vcf//g >  $downloadFile.sucess.temp
   grep -v -f $downloadFile.sucess.temp $downloadFile.temp > $output.temp
   echo -n "" >$output
   for line in `cat $output.temp`
   do
        grep $line $downloadFile >> $output
        grep $line $downloadFile
   done
   rm $output.temp
   rm $downloadFile.sucess*temp
}
#------------exe
# annovarAll 
# annovarVCF TCGA-BH-A0HX-01A-21D-A060-02.bam.var.vcf.gatk.vcf 
# filterAll 
# checkVcf2Download
# tempVCF="TCGA-E2-A14P-01A-31D-A19H-09.bam.var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf"
filtDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered/wu_rescue
tempVCF="TCGA-A7-A0CE-01A-11D-A12L-09.bam.var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf"
/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python-2.7.5/python $srcDir/filterMutation_local.py -i $tempVCF -o $filtDir/$tempVCF.filtered.vcf &  
# tempVCF="TCGA-A2-A04T-01A-21D-A128-09.bam.var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf"
# /ifs/home/c2b2/ac_lab/jh3283/tools/python/Python-2.7.5/python $srcDir/filterMutation_local.py -i $tempVCF -o $filtDir/$tempVCF.filtered.vcf 
