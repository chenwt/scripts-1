#!/bin/bash
#! -cwd
#By: J.He
#Desp.: combine gene list and miscellic operateions
#TODO: 

# awk "NR==1 || FNR !=1" brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt brca_geneSamplelist_UCceRNET_cancerGene_CNVMethFree_02072014.txt |sortxh|uniq  > brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt
sort brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt |uniq >brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt.temp
mv brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt.temp brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt
