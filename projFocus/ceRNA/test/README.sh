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

#head ~/SCRATCH/projFocus/ceRNA/data/snpArray/brca_snp_level3_839.mat > test_snp.mat
#~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/annot_SNP.py -i test_snp.mat -d ~/SCRATCH/database/projFocusRef/annot_snpArray_6.0_txt.bed  -o test_snp.mat.anno 

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


head -20 /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projFocus/ceRNA/data/methy/brca_meth_matrix_l3.mat.anno > input_test_meth_all.txt
~/tools/python/Python_current/python do_Beta2M_meth.py input_test_meth_all.txt 
mv do_Beta2M_meth.py ../do_Beta2M_meth.py
