#!/bin/bash
#By: J.He
#TODO: 



#/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python linkSNP_EXP.py
#/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python linkSNP_EXP.py -s test.snp -e test.exp -o test.out
#cp linkSNP_EXP.py uTest_snp_exp.py

#touch junk1 junk2 junk3
#./renameFiles.sh test1

# test for makeMatrnaseq.py level3
#i=0
#for file in `ls ~/SCRATCH/projFocus/ceRNA/data/rnaSeq/TCGA* |head -5`
#do
#  let i++
#  head -20 $file > test_$i.txt 
#  echo "test_$i.txt" >> input_test.txt
#done
#

#~/tools/python/Python_current/python makeMatrnaseq.py 
#~/tools/python/Python_current/python makeMatrnaseq.py -i input_test.txt -o output_test.txt 
#mv makeMatrnaseq.py ../makeMatRnaseql3.py


#mv ../get_RNAseqExpMat.py makeMat_General.py

#mv makeMat_General.py ../makeMat_General.py
#cp ../makeMat_General.py makeMat_General.py

#~/tools/python/Python_current/python makeMat_General.py -i input_test.txt -c 2 -e 2 -o output_test.txt 

#rm input_test.txt
#for file in `ls ~/SCRATCH/projFocus/ceRNA/data/snpArray/TCGA* |head -5`
#do
#  let i++
#  head -20 $file > test_$i.txt 
#  echo "test_$i.txt" >> input_test.txt
#done
#~/tools/python/Python_current/python makeMat_SNP.py -i input_test.txt -o output_test
#mv makeMat_SNP.py ../
#readlink -f ../makeMat_SNP.py
#rm *test*
####-------------test_for_anno_SNP.py-----------
#head ~/SCRATCH/projFocus/ceRNA/data/snpArray/brca_snp_level3_839.mat > test_snp.mat
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/annot_SNP.py -i test_snp.mat -d ~/SCRATCH/database/projFocusRef/annot_GenomeWideSNP_6_5cols_clean.txt -o test_snp.mat.anno 

#head -5 ~/SCRATCH/database/projFocusRef/annot_snpArray_6.0_txt.bed >test_input.txt
#sed -i 1d test_input.txt
#sed -i "1iID\tChr\tPos\tStrand" test_input.txt
#cat test_input.txt	

#awk '$2==22{print $0}' ~/SCRATCH/projFocus/ceRNA/result/snp/brca_snp_level3_839.mat.anno > input_test_snp.txt
#head -2 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_snp_level3_839.mat.anno > input_test_snp.txt
#head -2 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/edgeR_brca_l3DEGMat.txt.anno > input_test_exp.txt 

#~/tools/python/Python_current/python compareFileCols.py

#head -1 output_test_exp.txt |tr "\t" "\n" |wc -l 
#head -1 output_test_exp.txt |tr "\t" "\n" |head -10 
#head -1 output_test_snp.txt |tr "\t" "\n" |wc -l
#head -1 output_test_snp.txt |tr "\t" "\n" |head -10
#~/tools/python/Python_current/python compareFileCols.py -s input_test_snp.txt -e input_test_exp.txt
#rm *test*

#head -20 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/methy/brca_meth27k_matrix_l3.mat > input_test_meth27.txt
#head -20 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/methy/brca_meth450k_matrix_l3.mat > input_test_meth450.txt
##~/tools/python/Python_current/python compareFileRows.py  
#~/tools/python/Python_current/python compareFileRows.py -i input_test_meth27.txt -m input_test_meth450.txt -o output_test_meth.mat 
#cat output_test_meth.mat.log
#wc -l input_test_meth27.txt
#awk 'NR==1{print NF}' input_test_meth27.txt
#wc -l input_test_meth450.txt
#awk 'NR==1{print NF}' input_test_meth450.txt
#wc -l output_test_meth.mat
#awk 'NR==1{print NF}' output_test_meth.mat
#mv compareFileRows.py ../
#rm *test*

#head -20 /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projFocus/ceRNA/data/methy/brca_meth_matrix_l3.mat.anno > input_test_meth_all.txt
#~/tools/python/Python_current/python do_Beta2M_meth.py input_test_meth_all.txt 
#mv do_Beta2M_meth.py ../do_Beta2M_meth.py
#

#~/tools/python/Python_current/python make_Mat_cnv_level3.py -i input_test_cnv_tu.txt -o output_test.txt 
#mv make_Mat_cnv_level3.py ../


#####------test for make_Mat_cnvByGene_level3.py------
####
#head -20 /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projFocus/ceRNA/result/exp/edgeR_brca_l3DEGMat.txt.anno > test_exp.mat.anno
getNcol2Test() {
  file=$1
  numCol=9
  #awk -F"\t" 'BEGIN{OFS="\t"}{ out=$1;for(i=2; i<=9; i++){out=out"\t"$i};print out }' $file > $file.9cols
  awk -F"\t" 'BEGIN{ORS="\t"}{for(i=1; i<=9; i++) print $i;print "\n" }' $file > $file.9cols

}
##getNcol2Test test_meth_chr22.mat
#getNcol2Test test_exp.mat.anno 
##getNcol2Test test_snp_chr22.mat 
##sed -i "1isnpID\tchrom\tpos\tstrand\tsample1\tsample2\tsample3\tsample4\tsample5" test_snp_chr22.mat.9cols 
##sed -i "1imethID\tchrom\tpos\tstrand\tsample1\tsample2\tsample3\tsample4\tsample5" test_meth_chr22.mat.9cols 
#sed -i "1igeneid\tchrom\tpos\tstrand\tsample1\tsample2\tsample3\tsample4\tsample5" test_exp_chr22.mat.9cols 
#
##~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/merge_snp_meth_cnvl3_v1.py -s test_snp_chr22.mat.9cols -m test_meth_chr22.mat.9cols -c input_test_cnv_tu.txt -g test_exp_chr22.mat.9cols -o test_output_Gene_snp_meth_cnv.mat
#
##------------test_for_get_gene_cnv_mat------------
#head -5 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/CNV_snparray/level3/input_getMat_tu.txt > input_test_cnv_tu.txt
#cnt=0
#for line in `cat input_test_cnv_tu.txt`
#do
#  let cnt++
#  cat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/CNV_snparray/level3/$line > test_cnv_$cnt.txt 
#done
##~/tools/python/Python_current/python make_Mat_cnv_level3.py 
#ls test_cnv_?.txt |tr "\t" "\n" > input_test_cnv_tu.txt
#
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/mergeExpCNV.py -c input_test_cnv_tu.txt -e test_exp.mat.anno.9cols -o test_output_Gene_cnv.mat

#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/make_Mat_cnvByGene_level3.py -h  
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test/make_Mat_cnvByGene_level3.py -h  
#mv merge_snp_meth_cnvl3_v1.py ../
#------------------------------

####-------------test_for_prepareing_files-------
#~/tools/python/Python_current/python getSingleTSS_exp.py test_exp_chr22.mat.9cols test_exp_chr22.mat.9cols.uniTSS 
#~/tools/python/Python_current/python getCols.py -i test_exp_chr22.mat.9cols  -c input_test.txt -o output_test_exp_chr22_newcols.txt

###------test_for_grprep.r-------
#testinput is test_output_Gene_snp_meth_cnv.mat
#mv test_output_Gene_snp_meth_cnv.mat input_test_gene_snp_meth_cnv.mat
#mv test_exp_chr22.mat.9cols input_test_exp_chr22.mat.9cols

###----test_for_uTest----with_fdr_correction---
#head -10 /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno  > test.exp.anno 
##getNcol2Test test.exp.anno
#head -10 /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projFocus/ceRNA/data/snpArray/brca_snp_tumor_731.mat.anno > test.snp.anno 
##getNcol2Test test.snp.anno
#awk 'BEGIN{OFS="\t"}{print $1,"1",$3,$4,$5,$6,$7,$8,$9}' test.exp.anno > test.exp.an
#awk 'BEGIN{OFS="\t"}{print $1,"1",$3,$4,$5,$6,$7,$8,$9}' test.snp.anno > test.snp.an

#~/tools/python/Python_current/python filterSNP_utest_KWtest.py -e test.exp.an -s test.snp.an -o test_out_utest_snp_gene -j 1e-8
#~/tools/python/Python_current/python filterSNP_utest_KWtest_v2.py -e test.exp.an -s test.snp.an -o test_out_utest_snp_gene -j 0.2
#~/tools/python/Python_current/python filterSNP_utest_KWtest_v3.py -e test.exp.an -s test.snp.an -o test_out_utest_snp_gene -j 0.01
#head -100 /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/linkSNP_EXP/out_uTest_brca_snp_exp_pair.txt > test.input.pvalAdj.txt
#Rscript pvalAdj.r test.input.pvalAdj.txt test.out.pvalAdj.txt


####-------------------test for regression
#awk 'NR==1||$1=="ESR1"{print $0}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp/brca_exp_l3_731_DEG.mat.singleTSS.anno > input_test_reg_exp.mat
#awk 'NR==1||$1=="ESR1"{print $0}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/brca_gene_DEG_cnv_731.mat > input_test_reg_cnv.mat
#awk 'NR==1||$1=="ESR1"{print $0}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr9_KWtest.mat.anno.adjPass_1e-06.mat > input_test_reg_snp.mat 
#awk '$1=="ESR1"{print $0}' /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/snp/brca_gene_snp_chr6_KWtest.mat.anno.adjPass_1e-06.mat >> input_test_reg_snp.mat 
##creat a fake snp mat data
#cp /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test_knowBRCA.gene_snp_meth_cnv test_input_grpLassoSNP.txt
#Rscript grpLassoSNP.r 

##do this test in scratch----
cp input_test_reg*.mat /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test/ 
#echo "/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP_test.r input_test_reg_snp.mat input_test_reg_exp.mat input_test_reg_cnv.mat 0 " | qsub -l mem=4g,time=:30: -N test_grpreg -cwd 
