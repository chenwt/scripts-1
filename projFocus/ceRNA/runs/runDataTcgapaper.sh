#!/bin/bash
#$ -cwd
#By: J.He
#TODO: 

##--------get cancer gene from TCGA publication
getGeneListFromPaper(){
  cd ~/SCRATCH/projFocus/ceRNA/data/tcgaPaper/
  for file in `ls *txt`
  do 
    /ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/getGeneName.py -i ${file} -o ${file}_temp.genes
     sort ${file}_temp.genes|uniq > ${file}.genes
     echo "$file -----done"
  done
  rm *temp*
  cat *genes|sort|uniq -c |sort -k 1nr  > tcga.16papers.gene
}

##------get cancer gene regulator in ceRNET
getCeRNETRegulator(){
  cd ~/SCRATCH/projFocus/ceRNA/data/tcgaPaper/
  awk '$1>1{print $2}' tcga.16papers.gene > tcga.16papers.gene.freqBg1
  ~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/step1-1_extractCenetRegulator.py -i tcga.16papers.gene.freqBg1 -d /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/interactom/brca_ceRNA_network.txt -o tcga.16papers.gene.freqBg1.brcaCeRNETRegulator >> step1-1_extractCenetRegulator.py.log
}

##------check CrRegCanGene differential expression
# download data from TCGA do edgeR DEG analysis
# cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/topdown 
# tar -xvzf brca_CNVDNAmethRNAseqSom_level3_02042014.tar.gz
##move data to /ifs/data/c2b2/ac_lab/jh3283/projFocus/02042014/

PYTHON=~/tools/python/Python_current/python
###--------preprocesing data
getExpMatfromRnaseqLevel3(){
    expDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/rnaseqlevel3
    temp1="input_renameFiles.txt.temp"
    cd $expDir
    # echo -n "" > $temp1
    # for line in `ls *gene.quantification.txt`
    # do
    #   echo $line|awk 'BEGIN {FS="[. ]";OFS="\t"};{print $0,substr($2,0,19)".txt"}' >> $temp1
    # done
    # ~/scripts/projFocus/ceRNA/renameFiles.sh $temp1 
    temp2="input_makeMat.txt.temp"
    #mv renameFiles.sh.log $expDir/$temp2 
    cp /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/$temp2 . 
    # rm $temp1
    # cd $expDir
    output=$1
    output=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/$output 
    $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/makeMatRnaseql3.py -i $temp2 -o $output > $output.log
  }
# getExpMatfromRnaseqLevel3 brca_exp_level3_02042014.mat
getExpMatNormal(){
    cd /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/rnaSeq
    ls *-11A-* |tr "\t" "\n" > input.getExpMat.temp
    output=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/brca_expNormal_level3_02042014.mat 
    $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/makeMatRnaseql3.py -i input.getExpMat.temp -o $output > $output.log
}
# getExpMatNormal 

###---filter samples with cnv for each gene
cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/cnvlevel2/test
getCNVMat(){
    annoGeneStartEnd=/ifs/scratch/c2b2/ac_lab/jh3283/database/projFocusRef/refset_gene_start_end.tsv.sorted.uniq_resortedCol_geneSingleStartEnd
    # cnvDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/cnvlevel2
    cnvDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/02042014/cnvlevel2
    cd $cnvDir
    genelist=$1
    output=$2
    while read gene
    do 
      awk -v g=$gene '$1==g{print $0}' $annoGeneStartEnd >> genelist.geneSingleStartEnd
    done < $genelist
    # temp1="input_renameFiles.txt.temp"
    # grep "\.hg19.seg.txt" FILE_SAMPLE_MAP.txt |awk '{print $1"\t"$2}' > $temp1
    # ~/scripts/projFocus/ceRNA/step1-2.1_softlinkFiles.sh $temp1  
    $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-3_getGeneCNVMatLevel2.py -f step1-2.1_softlinkFiles.sh.log -g genelist.geneSingleStartEnd -o $output.temp 
    echo "transposing file..."
    if [ -f $output.temp ] ;then
      ~/bin/trfile $output.temp $output
    fi
    echo "cleaning temp files..."
    # rm step1-2.1_softlinkfiles.sh.log genelist.genesinglestartend input_renamefiles.txt.temp
    # ~/bin/rmlns $cnvDir
}
genelist='/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/tcgaPaper/tcga.16papers.genename' #the candidate cancer genes
output=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level2_02062014.mat
#getCNVMat $genelist $output

getCNVmatOld(){
    cnvDir=/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/CNV_snparray/level3
    cd $cnvDir
    output=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/cnv/brca_cnvTumor_level3_11112013.mat
    genelist='/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cancer.genelist.startend_022014.tsv'
    $PYTHON /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/step1-3_getGeneCNVMatLevel2.py -f /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/CNV_snparray/level3/input_getMat_tu.txt  -g ${genelist}  -o $output.temp
     ~/bin/trfile $output.temp $output
 }

countSample(){
awk 'NR>1{cnt=0;for(i=1;i<=NF;i++){
	  if($i==0) cnt = cnt+1};
	  print $1,cnt}' $output > $output.geneSampleCount
	}

getCNVmatOld $output
countSample $output
##-----count sample
###----filter samples with meth for each gene

genMethMat(){
   methDir=
  
}
genelist='/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/cancer.genelist.startend_022014.tsv'
##run DEG analysis
# cd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/


echo "##---END----"
