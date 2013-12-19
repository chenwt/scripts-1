#!/usr/bin/Rscript
#J.HE
#input: gene_exp; gene_cnv; gene_som; gene_snp;
#output:
#Descripiton: This is to use ridge linear fitting mode to fitting cnv+somatic mutation, and cnv+somatic + snp,
# 			  calculating the significance of SNPs' contribution for gene expression
# 


getArgs = function(){
	#base function used to get commandline argument specified by -- or by sequence
	args = commandArgs(trailingOnly = TRUE)
	hh <- paste(unlist(args),collapse=' ')
	if(length(grep("--",hh)) > 0){
	listOptions <- unlist(strsplit(hh,'--'))[-1]
	optionsArgs <- sapply(listOptions,function(x){
	         unlist(strsplit(x, ' '))[-1]
	        })
	optionsNames <- sapply(listOptions,function(x){
	  option <-  unlist(strsplit(x, ' '))[1]
	})
	names(optionsArgs) <- unlist(optionsNames)
 }else{
 	optionsArgs = args
 }
	return(optionsArgs)
}

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes//ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr21/test")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  # setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr21/test")
}

args = getArgs()
usage = "-Usage: Rscript ridgeFtest.r --exp exp.mat --snp snp.mat --cnv cnv.mat --som somaticMutation.mat --plotFlag 1/0[optional:default 0 ]
		Example: "
if(length(args) < 4 || is.null(args)){
	print(usage)
	stop("Input parameter error!")
}
print(args)


setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args$snp
inputexp    = args$exp
inputcnv    = args$cnv
inputsom    = args$som
plotFlag    = args$plotFlag
outputcoeff = paste(cwd,"/grplasso_coeff_",sep="")
outputFtest = paste(cwd,"/grplasso_Ftest_",sep="")
print(paste("current working directory",getwd()))


###------------test-start-------
inputsnp    = "input_snp_COL6A2"
inputexp    = "input_exp_COL6A2"
inputcnv    =  "input_cnv_COL6A2"
inputsom    = args$som
plotFlag    = args$plotFlag
outputcoeff = paste(cwd,"/grplasso_coeff_",sep="")
outputFtest = paste(cwd,"/grplasso_Ftest_",sep="")
print(paste("current working directory",getwd()))

###------------test-end--------


# load data
allData = getData(inputexp=inputexp,inputsnp=inputsnp,inputcnv=inputcnv)
dataExp = allData$dataExp
dataSnp = allData$dataSnp
