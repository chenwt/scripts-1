#!/usr/bin/Rscript
#J.HE
#input: gene_exp; gene_cnv; gene_som; gene_snp;
#output: file of gene snp contribution to expression significance
#Descripiton: This is to use ridge linear fitting mode to fitting cnv+somatic mutation, and cnv+somatic + snp,
# 			  calculating the significance of SNPs' contribution for gene expression
# Path: /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/test


sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  #setwd("/Volumes//ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr21/test")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("loaded from Linux")
  #setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr21/test")
}

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

args = getArgs()
usage = "-Usage: Rscript ftest.r --exp exp.mat --snp snp.mat --cnv cnv.mat --som somaticMutation.mat --plotFlag 1/0[optional:default 0 ]
		Example: "
if(length(args) < 5 || is.null(args)){
	print(usage)
	stop("Input parameter error!")
}

setwd(system("pwd",intern=T))
cwd         = getwd()
gene        = args['gene']
inputsnp    = args['snp']
inputexp    = args['exp']
inputcnv    = args['cnv']
inputsom    = args['som']
inputtype   = args['type'] #1,2,3
output      = paste(cwd,"/fTest_pval.txt",sep="")
print(paste("current working directory",getwd()))


# ###------------test-start-------
# cwd         = getwd()
# inputsnp    = "input_snp_COL6A2"
# inputexp    = "input_exp_COL6A2"
# inputcnv    =  "input_cnv_COL6A2"
# inputsom    =  "input_som_COL6A2"
# plotFlag    = 1
# gene        = "COL6A2"
# output = paste(cwd,"/fTest_pval.txt",sep="")
# print(paste("current working directory",getwd()))
# ###------------test-end--------

# load data
allData = getAllData(inputexp=inputexp,inputsnp=inputsnp,inputcnv=inputcnv,inputsom=inputsom)
allMergeData = allData$data_merge
require(MASS)
#-----model-exp=snp+cnv+som
calFtest = function(data,type="snp",outputFile=output,gene=gene){
  ##type: som, cnv, indel, snp
  switch(type,
         snp={
          	fit1 = lm(data[,1] ~ 1+ data[,'cnv'] + data[,'som'])    
            fit2 = lm(data[,1] ~ 1+data[,-1])
            ##anova sensitive to normal distribution	
            fresult =anova(fit1,fit2)
            names(fresult) = c("ResDf","RSS","DF","SumOfSquare","F","PrGreaterThanF")
         },
         som={},
         cnv={},
         indel={}
         )
  result = fresult$PrGreaterThanF[2]
  names(result) = paste(gene,type,sep="_")
  write.table(as.matrix(result),file=output,col.names=F,quote=F,sep="\t",append=T)
}
calFtest(allMergeData, outputFile=output, gene=gene)
#---model-exp=snp+cnv

#---model-exp=snp+som

#--add indel in the future

