#!/bin/bash
#$ -l mem=8g,time=8:: -cwd
#By: J.He
#TODO: 


####-------------function------

combineLog(){
 out=$1
 echo -e "chr\tnumGene\tnumTest\tnumPass" > $out
 for f in `ls *mat.anno.log`
 do
   chr=`echo $f|awk 'BEGIN{FS="_"}{gsub("chr","",$4);print $4}'`
   cntG=`awk 'BEGIN{FS=OFS="\t"}NR==5{print $2}' $f`
   cntpair=`awk 'BEGIN{FS=OFS="\t"}NR==6{print $2}' $f`
   cntPass=`awk 'BEGIN{FS=OFS="\t"}NR==7{print $2}' $f` 
   echo -e "$chr\t$cntG\t$cntpair\t$cntPass" >>$out
 done
 echo "#--END---"
}

combineSnpMatrix4Genes() {
 #head -1 brca_gene_snp_chr1_KWtest.mat.anno.adjPass_0.01.mat_GWASgene.mat > brca_gene_snp_KWtest.mat.anno.adjPass_0.01.mat_GWASgene.mat
 head -1 brca_gene_snp_chr1_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat > brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat
 for file in `ls *chr*GWASgene.mat`
 do
   awk 'NR>1{print $0}' ${file} >> brca_gene_snp_KWtest.mat.anno.adjPass_1e-06.mat_GWASgene.mat 
   echo -e "#${file}-----done---"
 done
 rm *chr*GWASgene.mat
}


##--------------function-edn
#ln -s ../../data/snpArray/brca_snp_level3_839.mat brca_snp_level3_839.mat
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/annot_SNP.py -i brca_snp_level3_839.mat  -d ~/SCRATCH/database/projFocusRef/annot_snpArray_6.0_txt.bed  -o brca_snp_level3_839.mat.anno

#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterSNP_utest_KWtest.py -s brca_snp_tumor_731.mat.anno -e /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno  -o  brca_gene_snp_KWtest.mat.anno 

#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v2.py -s brca_snp_tumor_731.mat.anno -e /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno  -o  brca_gene_snp_KWtest.mat.anno -j 1e-6 

#echo "~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731.mat.anno -e brca_exp_l3_731_DEG.mat.singleTSS.anno  -o  brca_gene_snp_KWtest.mat.anno -j 1e-6" |qsub -l mem=20g,time=48:: -cwd -N kwtestV3  

#python ~/HOME/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v2.py -s brca_snp_tumor_731.mat.anno -e brca_exp_l3_731_DEG.mat.singleTSS.anno  -o  brca_gene_snp_KWtest.mat.anno -j 1e-6 

###----get test data for snp_exp_KWTest
#awk '$2=="22"||NR==1{print $0}' brca_snp_tumor_731.mat.anno > brca_snp_tumor_731_chr22.mat.anno 
#awk '$2=="22"||NR==1{print $0}' brca_exp_l3_731_DEG.mat.singleTSS.anno > brca_exp_l3_731_DEG_chr22.mat.singleTSS.anno  
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731_chr22.mat.anno -e brca_exp_l3_731_DEG_chr22.mat.singleTSS.anno  -o  brca_gene_snp_chr22_KWtest.mat.anno -j 1e-6 

function chrFilterSNP(){
  chr=$1
  awk -v r=$chr '$2==r||NR==1{print $0}' brca_snp_tumor_731.mat.anno > brca_snp_tumor_731_chr$chr.mat.anno 
  awk -v r=$chr '$2==r||NR==1{print $0}' brca_exp_l3_731_DEG.mat.singleTSS.anno > brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno  
  echo "~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/test/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731_chr$chr.mat.anno -e brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno  -o  brca_gene_snp_chr${chr}_KWtest.mat.anno -j 1e-6" |qsub -l mem=8g,time=6:: -N kwtest${chr} -cwd  -e log/ -o log/ 
}

#chrFilterSNP "9"
#~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731_chr9.mat.anno -e brca_exp_l3_731_DEG_chr9.mat.singleTSS.anno  -o  brca_gene_snp_chr9_KWtest.mat.anno -j 1e-6
#chrFilterSNP "6"
#chrFilterSNP "1"
#chrFilterSNP "2"
#chrFilterSNP "3"
#chrFilterSNP "5"
#chrFilterSNP "7"
#chrFilterSNP "8"
#chrFilterSNP "11"
#chrFilterSNP "12"
#chrFilterSNP "13"
#chrFilterSNP "14"
#chrFilterSNP "15"
#chrFilterSNP "16"
#chrFilterSNP "17"
#chrFilterSNP "18"
#chrFilterSNP "19"
#chrFilterSNP "20"
#chrFilterSNP "21"
#chrFilterSNP "X"
#chrFilterSNP "Y"
#chrFilterSNP "10"

##need to fix the bug in python_snpannot
##for snpfile in `ls brca_gene_snp_chr*_KWtest.mat.anno.adjPass_1e-06.mat `
##do
##  echo "processing file: "$snpfile
##  sed -ic '1s/dbsnpID\t//g' $snpfile
##done


####-------------------------New code to run from Jan26,2014---
#get know snp annotated matrix
#grep -wf ../../knowledgeBase/GWAS_catalog_brca_SNPid.txt brca_snp_tumor_731.mat.anno > brca_snp_tumor_731_GWAS_DEG.mat.anno & 

function chrFilterSNP(){
   chr=$1
   cut=$2
   if [ ! -f brca_snp_tumor_731_chr$chr.mat.anno ]; then
	awk -v r=$chr '$2==r||NR==1{print $0}' brca_snp_tumor_731.mat.anno > brca_snp_tumor_731_chr$chr.mat.anno
   fi

   if [ ! -f  brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno ] ; then
     awk -v r=$chr '$2==r||NR==1{print $0}' brca_exp_l3_731_DEG.mat.singleTSS.anno > brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno
   fi
   #echo "~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/test/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731_chr$chr.mat.anno -e brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno -o brca_gene_snp_chr${chr}_KWtest.mat.anno -j 1e-2" |qsub -l mem=8g,time=6:: -N kwtest${chr} -cwd  -e log/ -o log/
   ~/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/test/filterSNP_utest_KWtest_v3.py -s brca_snp_tumor_731_chr$chr.mat.anno -e brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno -o brca_gene_snp_chr${chr}_KWtest.mat.anno -j $cut  &
      }


#for snpfile in `ls brca_snp_tumor_731_chr*.mat.anno`
#do 
#  head -1 brca_snp_tumor_731.mat.anno > test/${snpfile}_GWASCataSNP.anno
#  grep -wf ../../knowledgeBase/GWAS_catalog_brca_SNPid.txt brca_snp_tumor_731.mat.anno >> test/${snpfile}_GWASCataSNP.anno
#  echo -e "done with $snpfile"
#done 
#gene geneName for snps
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterSNP_utest_KWtest_v3.py -s /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/snpArray/brca_snp_tumor_731.mat.anno  -e /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/brca_exp_l3_731_DEG.mat.singleTSS.anno  -o brca_snp_tumor_KWtest.mat.anno -j 1e-2
# ###-----redo KW test on new snp matrix 

function chrFilterSNP_v2(){
 chr=$1
 cut=$2
 if [ ! -f brca_snp_tumor_731_chr$chr.mat.anno ] ; then
   echo -e "no snp matrix... please run model 1"
 fi
 if [ ! -f brca_exp_l3_731_DEG_chr$chr.mat.singleTSS.anno ] ; then
   echo -e "no snp matrix... please run model 2"
 fi
 awk '{}' 
 ~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/extractSNPmat.py -i temp$chr.genelist -m brca_snp_tumor_731_chr$chr.mat.anno & 
}

#chrFilterSNP '1' 1e-6 &
# for chr in  `seq 7 22` Y
# do
#  chrFilterSNP "$chr" 1e-6
# done 

###---------extracting some gene's SNP matrix
# for file in `ls brca_gene_snp_chr*_KWtest.mat.anno.adjPass_0.01.mat`
 # for file in `ls brca_gene_snp_chr*_KWtest.mat.anno.adjPass_1e-06.mat`
 # do
 #     ~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/extractSNPmat.py -i /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_allGeneName.txt -m /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/$file -o ${file}_GWASgene.mat
 #     echo -e "${file} done...."
 # done 

##====combining files

# combineSnpMatrix4Genes
# combineLog "kwtestSummary.txt"
wc -l brca_gene_snp_chr*KWtest.mat.anno.adjPass_0.01.mat >> kwtestSummary_0.01.txt
