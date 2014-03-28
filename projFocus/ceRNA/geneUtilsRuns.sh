#!/bin/bash
#! -cwd
#By: J.He
#Desp.: all functions for RUNs in projFocus, general usage
#TODO: 

CDT=`date|awk '{print $2"-"$3"-"$6}'`
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA
PYTHON=~/tools/python/Python_current/python
RSCRIPT=~/tools/R/R_current/bin/Rscript
###-----general
getGeneAnno(){
    ##annotate gene with gene transcription starting site and last exon end position
    annoGeneStartEnd=/ifs/data/c2b2/ac_lab/jh3283/database/refseq/refseq_gene_hg19_selected_Mar22_Tsstse.tsv.single.tsv
    fname=$1
    out=$2
    echo -n "" > $out 
    while read gene
    do 
      awk -v g=$gene '$1==g{print $0}' $annoGeneStartEnd >> $out 
    done < $fname
} 

###get Expression matrix
getExpMatfromRnaseqLevel3(){
    ### given the raw level3 file folder, output a matrix of raw counts  
    ### example getExpMatfromRnaseqLevel3 
    # expDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/rnaseqlevel3
    expDir=$1
    output=$2
    temp1=$expDir/"input_renameFiles.txt.temp"
    grep "gene.quantification.txt" $expDir/FILE_SAMPLE_MAP.txt |awk -v d=$expDir '{print d"/"$1"\t"d"/"substr($2,6,11)}' > $temp1
    ~/scripts/projFocus/ceRNA/step1-2.1_softlinkFiles.sh $temp1
    temp2=$temp1.log 
    cd $expDir
    $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/makeMatRnaseql3.py -i $temp2 -o $output > $output.log
    rm $temp1 $temp2
    ~/bin/rmlns $expDir
  }

getCNVMat(){
    ##---generate CNV matrix for each gene inputed, output CNV value for each gene
    ##input: cnv directory; genelist; output matrix name
    cnvDir=$1
    genelist=$2
    output=$3
    temp0=$genelist.temp
    getGeneAnno $genelist $temp0
    temp1=$cnvDir/"input_renameFiles.txt.temp"
    grep "\.hg19.seg.txt" $cnvDir/FILE_SAMPLE_MAP.txt |awk -v d=$cnvDir '{print d"/"$1"\t"d"/"substr($2,6,11)}' > $temp1 
     ~/scripts/projFocus/ceRNA/step1-2.1_softlinkFiles.sh $temp1  
    temp2=$nput_renameFiles.txt.temptemp1.log
    temp3=$output.temp
    $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-3_getGeneCNVMatLevel2.py -f $temp2 -g $temp0 -o $temp3 
    echo "transposing file..."
    if [ -f $temp3 ] ;then
      ~/bin/trfile $temp3 $output
    fi
    echo "cleaning temp files..."
    rm $temp0 $temp1 $temp2 $temp3 
    ~/bin/rmlns $cnvDir
}

countSample(){
      awk 'NR>1{cnt=0;for(i=1;i<=NF;i++){
	  if($i==0) cnt = cnt+1};
	  print $1,cnt}' $output > $output.geneSampleCount
	}

genMethMat(){
   ##changes at 02172014
   methDir=$1
   genelist=$2
   output=$3

   temp1=$methDir/"input_linkfile.txt.temp"
   grep "BRCA.HumanMethylation" $methDir/FILE_SAMPLE_MAP.txt |awk -v d=$methDir '{print d"/"$1"\t"d"/"substr($2,6,11)}' > $temp1 
   ~/scripts/projFocus/ceRNA/step1-2.1_softlinkFiles.sh $temp1 ## generate log file, which is used as input for makeMat script 
   temp2=$temp1.log 
   temp3=$output.temp
   $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-5_makeMethMatlevel3.py -f $temp2  -g $genelist -o $temp3 
    
   ~/bin/trfile $temp3 $output 
   # rm $temp0 $temp1 $temp2 $temp3 
    ~/bin/rmlns $methDir
}

genMethMat_v2(){
   ##changes at 02172014
   methDir=$1
   genelist=$2
   output=$3

   temp1=$methDir/"input_linkfile.txt.temp"
   grep "BRCA.HumanMethylation" $methDir/FILE_SAMPLE_MAP.txt |awk -v d=$methDir '{print d"/"$1"\t"d"/"substr($2,6,11)}' > $temp1 
   ~/scripts/projFocus/ceRNA/step1-2.1_softlinkFiles.sh $temp1 ## generate log file, which is used as input for makeMat script 
   temp2=$temp1.log 
   temp3=$output.temp
   $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-5_makeMethMatlevel3_v2.py -f $temp2  -g $genelist -o $temp3 
    
   ~/bin/trfile $temp3 $output 
   # rm $temp0 $temp1 $temp2 $temp3 
    ~/bin/rmlns $methDir
}

##-----------get Target Gene somatic mutation matrix for filtering sample
getCGSomMat(){
  vcfDir='/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/wgsVars/filtered'
  somDir='/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som'
  cd $somDir
  $PYTHON  $srcDir/varcall/step-6_getMAF.py -d $vcfDir -r gene -k snp -g $genelist -o $somDir/brca_somAll_targetGeneMut_2014-03-01.mat 
}


echo  "~/scripts/projFocus/ceRNA/geneUtilsRuns.sh loaded"
