#!/bin/bash
#! -cwd
#By: J.He
#Desp.
#TODO: 

###test one on ESR1
PYTHON=~/tools/python/Python_current/python
testDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/
# grep -w ESR1 /ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt|awk -F"\t" '{print $1"\n"$2}' |sort|uniq > ESR1.regulator
# grep -w -f ESR1.regulator /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulator.txt.tss > ESR1.regulator.deg.tss
# ~/tools/python/Python_current/python ~/bin/splitByKey.py -i /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat -s 500  -k /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/ESR1.regulator.deg.tss -o /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/brca_somTumorWU_ESR1_Regulator_v3_Mar-6-2014.mat_dist100_.mat-temp

input=$somDir/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat
genelist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulator.txt.tss
barcode=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/brca_somTumorWU_barcode.list
getMutMatByRegion=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getMutMatByRegion_ing.py

output=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/brca_somTumorWU_ESR1_Regulator_v3_Mar-6-2014.mat_dist100_.mat
dist=100
###splitting step has been done
# $PYTHON $getMutMatByRegion -i $input -g $genelist -s $barcode -d $dist -o $output  
# mv $output /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/brca_somTumorWU_ESR1_Regulator_v3_Mar-6-2014.mat_dist100.mat
# dist=1000
# $PYTHON $getMutMatByRegion -i $input -g $genelist -s $barcode -d $dist -o $output  
# mv $output /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/brca_somTumorWU_ESR1_Regulator_v3_Mar-6-2014.mat_dist1000.mat

dist=0
# $PYTHON $getMutMatByRegion -i $input -g $genelist -s $barcode -d $dist -o $output  
# mv $output /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test/brca_somTumorWU_ESR1_Regulator_v3_Mar-6-2014.mat_dist0.mat


#####----test on APC
gintDegSmplist=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt
crtDir=/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/test
srcDir=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getMutMatByRegion_ing.py

g="APC"

# grep -w $g $gintDegSmplist > $crtDir/APC.smps 
# grep -w $g /ifs/data/c2b2/ac_lab/jh3283/projFocus/other/brca_ceRNA_network.txt|awk     -F"\t" '{print $1"\n"$2}' |sort|uniq > $crtDir/$g.regulators
# grep -w -f $crtDir/$g.regulators /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulator.txt.tss > $crtDir/$g.regulator.deg.tss


# ~/tools/python/Python_current/python ~/bin/splitByKey.py -i /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat -s 500  -k $crtDir/$g.regulator.deg.tss -o $crtDir/brca_somTumorWU_${g}_Regulator_v3_Mar-6-2014.mat-temp

# splitDir=$crtDir/brca_somTumorWU_APC_Regulator_v3_Mar-6-2014.mat-temp/
# getMutMatByRegion=/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/getMutMatByRegion_ing.py
# dist=1000
# output=$crtDir/brca_somTumorWU_APC_Regulator_${cdt}_dist${dist}.matrix
# cdt=`date|awk '{print $2"-"$3"-"$6}'`
# ~/tools/python/Python_current/python  $getMutMatByRegion -i $splitDir -d $dist -s $barcode -o $output 


###----test---
# awk -F"\t" 'NR==1||$1=="KIF23"{print $0}' ../brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist0_Mar-7-2014.matrix > $crtDir/CHEK1_KIF23.som.mat

# awk '($3>=69866056&&$4<=69866405)||($3>=70052447&&$4<=70052509)' KIF23_mut_raw.txt >> KIF23_mut_keyPoint.txt 
# for key in `cat key_name.txt`; do grep $key KIF23_4colraw.txt | awk 'BEGIN{FS="\t|:"}$3>=69380810&&$3<69395000{print $0}' >>key_name.txt.mut; done 
