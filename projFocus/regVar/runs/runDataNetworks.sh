#!/bin/bash
#$ -cwd
#By: J.He
#Desp.:

##--get 4col regulon, singulon
Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/regVar/step0/step0_getTxtAracneRegulon.r

##--ceRNET to entrez id
Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/regVar/step0/step0_sym2entrezid_cernet.r



