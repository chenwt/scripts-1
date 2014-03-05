#!/bin/bash
#! -cwd
#By: J.He
#Desp.: combine gene list and miscellic operateions
#TODO: 

# awk "NR==1 || FNR !=1" brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt brca_geneSamplelist_UCceRNET_cancerGene_CNVMethFree_02072014.txt |sortxh|uniq  > brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt
# sort brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt |uniq >brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt.temp
# mv brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt.temp brca_geneSamplelist_nature12912NovelGene_CNVMethFree_02182014.txt.deg_02102014.txt_regulator.txt
# awk -F"\t" '{print $1}'  brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt.deg_2014-02-24.txt |sort |uniq > brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt.deg_2014-02-24.txt.genelist 
# ~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/annotGeneByStartEndPos.py -i brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt.deg_2014-02-24.txt.genelist -o brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt.deg_2014-02-24.txt.genelist.singleStartEnd 

~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/annotGeneByStartEndPos.py -i brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt_regulator.txt -o brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt_regulator.txt.singleStartEnd 
# awk -F"\t" '{n=split($2,a,";")
#     if (n>10)
#       print $0
#     }' brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt  > brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.10moreSmps


