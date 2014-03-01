#!/bin/bash
#! -cwd
#By: J.He
#Desp.: given CGgeneSample list, ceRnet, extract regulator and samples list
#TODO:

PYTHON=/ifs/home/c2b2/ac_lab/jh3283/tools/python/Python_current/python 
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/
gslist="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_geneSamplelist_combined_CG_CNVMethFree_02242014.txt.deg_2014-02-24.txt"
cernet="/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_ceRNA_network.txt"
out=${gslist}_regulatorSamples
# head -1 $gslit > $out
$PYTHON $srcDir/processData/getRegulatorSamplelist.py -i $gslist  -d $cernet  -o $out 

