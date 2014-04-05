#!/bin/bash
#By: J.He
#TODO: 

CWD=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/indel
PYTHON=~/tools/python/Python_current/python
##--making bamlist files--------------------------------------------
###----------------make--input-list-file--end--------------------------

#-------function-to-submit---bam-file-list---------
qsubRun (){
  ##----input: bam file list with full path
  cnt=0
  while read line
  do
    let cnt=$cnt+1
    jobid=`echo $line|awk 'BEGIN{FS="/"}{split($13,a,"-");print a[3]"-"a[4]}'`"_$cnt"
    # echo $jobid
    echo "/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/indelCall_v5.sh $outDir $line"|qsub -l mem=8g,time=160:: -N $jobid -cwd -o ./log -e ./log >>qsub.log
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

getInput(){
  pidlist=$1
  echo -n "" >$pidlist.inputBamlist
  while read pid
  do
    tempBamArray=( `ls ${dataDir}/${pid}*.bam ` )
    for i in ${tempBamArray[@]}
    do
      readlink -f $i >> $pidlist.inputBamlist
    done
  done < $pidlist
}

batchQsubRun(){
  for i in `seq 1 8`
  do
    qsubRun pid.txt.inputBamlist_$i 
    sleep 48hs 
  done
}


###-------function----end--------

#----data folder:
dataDir=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/
outDir=$(pwd)
# getInput $CWD/pid.txt
# ~/bin/splitByN $CWD/pid.txt.inputBamlist 6
# qsubRun pid.txt.inputBamlist_2
# qsubRun pid.txt.inputBamlist_3
# qsubRun pid.txt.inputBamlist_4
# qsubRun pid.txt.inputBamlist_5
# qsubRun pid.txt.inputBamlist_6
# qsubRun pid.txt.inputBamlist_7
# qsubRun pid.txt.inputBamlist_8

###---------run-----using--input---file---list

# qsubRun pid.txt.inputBamlist_part1
# qsubRun pid.txt.inputBamlist_part2 
##------g

# touch input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt
# for inputlist in input.bam.list.4 input.bam.list.5 input.bam.list.6
# do
#   grep -f ../../data/sampleInfo/brca_wgs_bam_summary_02042014.tsv_TumorSample $inputlist >> input.bam.list.brca_wxsInwgsNotDone.bam.tumor_02022014.txt
# done 

# ###------rename all files----
# while read pid
# do
#   normal=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/14A/{print $16}'`
#   normal=${normal}_dindel_ouput.variantCalls.VCF.filtered.VCF
#   tumor=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/0[93]A/{print $16}'`
#   tumor=${tumor}_dindel_ouput.variantCalls.VCF.filtered.VCF
#   relapse=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/04A/{print $16}'`
#   relapse=${relapse}_dindel_ouput.variantCalls.VCF.filtered.VCF
#   out=$pid.intersect
#   # echo -e $pid"\t"$out"\n"$normal"\n"$tumor"\n"$relapse
#   ~/tools/python/Python_current/python  ~/scripts/projAML/WXS/indel/intersectIndel4one.py -n $normal -t $tumor -r $relapse -o $out
# done < pid.txt

###---annotation snpEFF

# for file in `ls $CWD/*variantCalls.VCF.filtered.VCF`
# do
#   out=$file.snpeff.vcf 
#   java -Xmx4g -jar /ifs/home/c2b2/ac_lab/jh3283/tools/snpEff/snpEff.jar eff -c /ifs/home/c2b2/ac_lab/jh3283/tools/snpEff/snpEff.config -v GRCh37.69 $file > $out 
#   $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/parseSnpeffInfo_v3.py $out Y 
# done 

# ##-----second filtering 
# while read pid
# do
#   normal=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/14A/{print $16}'`
#   normal=${normal}_dindel_ouput.variantCalls.VCF.filtered.VCF.snpeff.vcf.parsed...var
#   tumor=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/0[93]A/{print $16}'`
#   tumor=${tumor}_dindel_ouput.variantCalls.VCF.filtered.VCF.snpeff.vcf.parsed...var
#   relapse=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/04A/{print $16}'`
#   relapse=${relapse}_dindel_ouput.variantCalls.VCF.filtered.VCF.snpeff.vcf.parsed...var
#   out=$pid.2ndfilter.snpeff.intersect
#   annofile=~/DATA/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv 
#   $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/filterIndel.py -n $normal -t $tumor -r $relapse -a $annofile -o $out >> $CWD/filtet2nd.snpeff.log
# done < pid.txt

##----find recurrence genes
# awk '{print $1":"$2}' *2ndfilter.snpeff.intersect.relapse| sort |uniq -c |sort -k 1 -r  


# #-----------ln file
# vcfDir=/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/indel/raw
# while read pid
# do
#   normal=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/14A/{print $16}'`
#   normal=${normal}_dindel_ouput.variantCalls.VCF
#   # ln -s $vcfDir/$normal $CWD/${pid}_normal.vcf 

#   tumor=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/0[93]A/{print $16}'`
#   tumor=${tumor}_dindel_ouput.variantCalls.VCF
#   # ln -s $vcfDir/$tumor $CWD/${pid}_tumor.vcf 

#   relapse=`grep -w $pid pid.txt.inputBamlist |awk -F"/" '$13~/04A/{print $16}'`
#   relapse=${relapse}_dindel_ouput.variantCalls.VCF
#   # ln -s $vcfDir/$relapse $CWD/${pid}_relapse.vcf 
  

#   ##-----
#   cosmic=/ifs/data/c2b2/ac_lab/jh3283/database/cosmic/CosmicMutantExportCensus_v68.tsv.locus
#   normal=$CWD/${pid}_normal.vcf
#   normalout=$normal.cosmic.vcf
#   tumor=$CWD/${pid}_tumor.vcf
#   tumorout=$tumor.cosmic.vcf
#   relapse=$CWD/${pid}_relapse.vcf
#   relapseout=$relapse.cosmic.vcf
#   # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/getCosmicIndel.py -i ${normal} -c $cosmic -o $normalout
#   # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/getCosmicIndel.py -i ${tumor} -c $cosmic -o $tumorout
#   # $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/getCosmicIndel.py -i ${relapse} -c $cosmic -o $relapseout
#   pidCosmic=${pid}.cosmic.vcf
#   awk 'NR==1||FNR!=1{print FILENAME,"\t", $0}'  $normalout $tumorout $relapseout > $pidCosmic
# done < pid.txt

##-----------find cosmic mutations
###----get cosmic indels (short indel)

echo "#-----END---"

function pidFilter(){
  pid=$1
  pattern=$2
  dir=$3
  normal=`grep -w $pid $CWD/pid.txt.inputBamlist |awk -F"/" '$13~/14A/{print $16}'`
  normal=$dir/${normal}${pattern}

  tumor=`grep -w $pid $CWD/pid.txt.inputBamlist |awk -F"/" '$13~/0[93]A/{print $16}'`
  tumor=$dir/${tumor}${pattern}

  relapse=`grep -w $pid $CWD/pid.txt.inputBamlist |awk -F"/" '$13~/04A/{print $16}'`
  relapse=$dir/${relapse}${pattern}  
  annofile=~/DATA/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv 
  out=${pid}.snpeff.2ndftilter_mar31
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/filterIndel_v4.py  -r ${relapse} -t ${tumor} -n ${normal} -a $annofile -o $out

}

pid="PASFEW"
ptn="_dindel_ouput.variantCalls.VCF.filtered.VCF.snpeff.vcf.parsed...var"
dir=/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/indel/filter2nd
# while read pid 
# do 
#   echo "#-----${pid}-----"
#   pidFilter $pid $ptn $dir
# done < pid.txt  
while read pid 
# pid="PASFEW"
do
  tu=${pid}.snpeff.2ndftilter_mar31.tumor
  re=${pid}.snpeff.2ndftilter_mar31.relapse
  cosmic=/ifs/data/c2b2/ac_lab/jh3283/database/cosmic/CosmicMutantExportCensus_v68.tsv.locus
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/getCosmicIndel.py -i ${re} -c $cosmic -o $re.cosmic
  $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/indel/getCosmicIndel.py -i ${tu} -c $cosmic -o $tu.cosmic
done < pid.txt
