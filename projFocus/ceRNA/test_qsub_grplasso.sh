#!/bin/bash
#By: J.He
#TODO: 
#$ -l mem=4g,time=:30:
#$ -cwd

/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP_test.r input_test_reg_snp.mat input_test_reg_exp.mat input_test_reg_cnv.mat 0


