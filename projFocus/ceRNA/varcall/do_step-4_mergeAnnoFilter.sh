#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/varcall
filtDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered

function annovarAll(){
  rawDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars
  cd $rawDir
  rawVCF_array=( `ls TCGA-*var.vcf.gatk.vcf` )
   echo ${tempDir_array[@]}
  for tempVCF in "${rawVCF_array[@]}"
  do
      cntFiltVCF=`ls $tempVCF.annovar.summary.genome_summary.csv.vcf 2>/dev/null |wc -l`
      if [ $cntFiltVCF  == 0 ]
      then
      	  # echo "annovar annotating... $tempVCF"
	  tempVCFp=`readlink -f $tempVCF`
    	  pid=`echo $tempVCF|awk -F"-" '{print $3}' `
	  echo "$srcDir/step-5_doAnnovarAll.sh $tempVCFp " | qsub -l mem=4g,time=8:: -N af_$pid -e ./log -o ./log -cwd >>qsub.log
	  tail -1 qsub.log
      fi
  done
}

function filterAll(){
  rawDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/rawVars
  filtDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered
  cd $rawDir
  rawVCF_array=( `ls TCGA-*var.vcf.gatk.vcf.annovar.summary.genome_summary.csv.vcf` )
  echo ${tempDir_array[@]}
  for tempVCF in "${rawVCF_array[@]}"
  do
      temp=`echo $tempVCF|sed s/.annovar.summary.genome_summary.csv.vcf//g`
      if [ -f $temp ]; then 
	echo " deleting $temp" 
	rm $temp 
      fi
      cntFiltVCF=`ls $filtDir/$tempVCF.annovar.summary.genome_summary.csv.vcf* 2>/dev/null |wc -l`
      if [ $cntFiltVCF  == 0 ]
      then
	  echo "filter... $tempVCF"
    	  pid=`echo $tempVCF|awk -F"-" '{print $3}' `
	  echo " /ifs/home/c2b2/ac_lab/jh3283/tools/python/Python-2.7.5/python $srcDir/filterMutation_local.py -i $tempVCF -o $filtDir/$tempVCF.filtered.vcf" | qsub -l mem=4g,time=8:: -N af_$pid -e ./log -o ./log -cwd
      fi
  done
}


#------------exe
annovarAll 
