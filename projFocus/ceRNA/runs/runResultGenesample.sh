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

# ~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/annotGeneByStartEndPos.py -i brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt_regulator.txt -o brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.deg_20140303.txt_regulator.txt.singleStartEnd 
# awk -F"\t" '{n=split($2,a,";")
#     if (n>10)
#       print $0
#     }' brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt  > brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt.10moreSmps

# awk -F"\t" 'NR>1{n=split($2,a,";"); for (i=1;i<=n;i++) print a[i];}' brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt |sort|uniq|less 

# echo -e "Target\tNumOfRegulator" > brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulatorFreq.txt

# while read line  
# do 
#   g=`echo $line|awk '{print $1}'`
#   gCnt=`awk -F"\t" -v g=$g '$1==g||$2==g{print $1"\n"$2}' /ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt |sort|uniq |awk -v g=$g '$1!=g' |wc -l 2>/dev/null` 
#   echo -e $g"\t"$gCnt >> brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulatorFreq.txt
# done < brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt 

# awk '$2>0{print $1}' brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulatorFreq.txt > temp
# grep -w -f temp brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt > brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt.hasRegulator.txt 
# rm temp
echo "#---DONE---"
