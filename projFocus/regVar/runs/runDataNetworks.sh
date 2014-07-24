#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

##--get 4col regulon, singulon
# Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/regVar/step0/step0_getTxtAracneRegulon.r

##--ceRNET to entrez id
# Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/regVar/step0/step0_sym2entrezid_cernet.r

tfnet=/ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_geneid-regulon.rda.4col.txt.naomit.symbol
~/tools/python/Python_current/python /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/processData/pickleDumpCernet.py $tfnet

