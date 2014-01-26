#!/bin/bash
#By: J.He
#TODO: 
#mkfifo pipein pipeout
##ln -s ../../data/rnaSeq/brca_rnaseq_rawCount_l3.mat brca_rnaseq_rawCount_l3.mat 
#ln -s ../../data/rnaSeq/brca_rnaseq_rawCount_l3_geneSymbol.mat brca_rnaseq_rawCount_l3_geneSymbol.mat 
##prepare for the design file
#grep -f ../../data/rnaSeq/input_makeMatRnaseql3.txt ../../data/rnaSeq/input_renameFiles.txt > pipein &
#cat pipein| awk 'BEGIN{OFS="\t"}{split($2,a,"-");print $2,$1,a[3],substr(a[4],0,2)}' > input_edgeR_ColDesign.txt
#
#cp ../grpreg/final.header_exp.sorted .
####-------Step2---getting-DEG-normalized-dataset---need R3-02--
#/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript ~/scripts/projFocus/ceRNA/getDEG_Mat_edgeR.r /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp input_edgeR_ColDesign.txt brca_rnaseq_rawCount_l3_geneSymbol.mat brca_rnaseq_rawCount_l3_geneSymbol final.header_exp.sorted 

###----step_3----annotate---genes
##/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/annot_SNP.py -i edgeR_brca_l3DEGMat.txt -d /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/database/projFocusRef/annot_entrezID.bed -o edgeR_brca_l3DEGMat.txt.anno 
#/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python ~/scripts/projFocus/ceRNA/annot_SNP.py -i brca_rnaseq_rawCount_l3_geneSymbol_DEG_Mat.txt -d#/ifs/scratch/c2b2/ac_lab/jh3283/database/projFocusRef/refseq_gene_tss_geneSymbol.signleTSS.5cols -o brca_exp_l3_731_DEG.mat.singleTSS.anno 
#rm pipein pipeout

#sed -ic '1s/strand\tGene/strand/g' brca_exp_l3_731_DEG.mat.singleTSS.annoc 
#cut -f1 brca_exp_l3_731_DEG.mat.singleTSS.anno |sed 1d > brca_exp_l3_731_DEG.mat.signleTSS.anno.GeneName



####-------------------------run after Jan.2014 for known snp/genes
while read line 
do
grep -w $line brca_exp_l3_731_DEG.mat.singleTSS.anno >> brca_exp_l3_731_DEG.mat.singleTSS.anno.GWASCataGene.mat.anno

done < /ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/knowledgeBase/GWAS_catalog_brca_GeneName.txt 
