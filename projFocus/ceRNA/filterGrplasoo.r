##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: output files prefix from grpLassoSNP.r > 
#output: <r data file: >
#Usage: Rscript grplassoSNP.r input
#Description: this file was created for projFocus, ceRNA, used after running grpLassoSNP.r to filter significant gene-snps 
#TODO: develop plot option!

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test")
}

###-----------------functions-----------------

###------------------functions---end---------

#getting command line parameters
args = getArgs()
usage = "Usage: Rscript filterGrplasso.r --exp exp.mat --snp snp.mat --cnv cnv.mat --som som.mat --out outputfile --type 1[2/3,model type] [optional: --plot 1/0[plot/no plot] --nperm 1000]"
example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r "
if(length(args) < 7 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop("Input parameter error!")
}

setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args['snp']
inputexp    = args['exp']
inputcnv    = args['cnv']
inputsom    = args['som']
type        = args['type']
plotflag    = as.integer(args['plot']) #optional 
# nperm       = as.integer(args['nperm'])
output      = paste(cwd,"/",args['out'], sep="")
print(paste("current working directory:",cwd))
# print(paste(inputsnp,inputexp,inputcnv,inputsom,type,plotflag,nperm,output,sep="\n"))
##-----------------------test-##----------------------------
# setwd(system("pwd",intern=T))
# cwd         = getwd()
# inputsnp    = "input_test_reg_snp.mat"
# inputexp    = "ESR1_brca_exp_l3_731_DEG.mat.singleTSS.anno"
# inputcnv    = "ESR1_brca_gene_DEG_cnv_731.mat"
# inputsom    = "ESR1_brca_somForDeg.mat"
# type        = 1
# plotflag    = 1
# nperm       = 100
# output      ="grplasso_coeff_"
# output      = paste(cwd,"/",output, sep="")
# print(paste("current working directory:",cwd))
##----------------------------##----------------------------