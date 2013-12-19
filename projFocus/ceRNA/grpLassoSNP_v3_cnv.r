##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: expression and snp gt for one gene > 
#output: <file: gene: snp with weight 
#Usage: Rscript grplassoSNP.r input
#Description: This is used in projFocus/ceRNA/grpreg, to select SNPs contributing to gene differential expression
#             require package irr, grpreg,gplots
#             Major change: with transformation, 
#                           functionlized all main parts, 
#                           residual calculation corrected
#TODO:add quality checking script for each file and the content

#getting command line parameters
plotflag    = 0 # "0  for no plotting "
args        = commandArgs(TRUE)
if (is.null(args) || length(args) < 3){
  print("Please provide parameters in right order!")
  print("Usage: Rscript grplassoSNP_v1.r <inputsnp> <inputexp> <inputcnv> <1/0>")
  stop("Stop since no parameters")
}else{
  print(args)
}

setwd(system("pwd",intern=T))
cwd         = getwd()
inputsnp    = args[1]
inputexp    = args[2]
inputcnv    = args[3]
plotflag    = args[4]
outputcoeff = paste(cwd,"/grplasso_coeff_",sep="")
outputFtest = paste(cwd,"/grplasso_Ftest_",sep="")
print(paste("current working directory",getwd()))
# print(paste("Current working directory:",system("pwd",intern=T)))


##-------load useful functions-------------------
sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoFunctions.R")
  # setwd("/Volumes//ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr21/test")
  #setwd("/Users/jh3283/projFocus") 
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoFunctions.R")
  # setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/chr21/test")
}

##----------------------------##----------------------------
##-----------------------test-----------------------------
# outputcoeff = "grplasso_coeff_"
# outputcoeff = "grplasso_Ftest_"
#  inputsnp    = "temp/input_snp_COL6A2"
#  inputexp    = "temp/input_exp_COL6A2"
#  inputcnv    = "temp/input_cnv_COL6A2"
#  system(paste("cp  ../",inputsnp,"  .",sep=""))
#  system(paste("cp  ../",inputexp,"  .",sep=""))
#  system(paste("cp  ../",inputcnv,"  .",sep=""))
#  system(paste("cp ../","grplasso_COL6A2.RData"," .",sep=""))
# inputsnp    = "input_snp_COL6A2"
# inputexp    = "input_exp_COL6A2"
# inputcnv    = "input_cnv_COL6A2"

# # inputsnp    = "input_test_reg_snp.mat"
# # inputexp    = "input_test_reg_exp.mat"
# # inputcnv    = "input_test_reg_cnv.mat"
##----------------------------##----------------------------


##----------------------------##----------------------------
getData = function(inputexp,inputsnp,inputcnv,plotflag){
  ##load exp data
  dataExp             = read.table(inputexp,header =T)
  rownames(dataExp)   = dataExp[,1]
  genename            = dataExp[,1] 
  dataExp             = sapply(dataExp[,-c(1:4)],
                              function(x){as.numeric(as.character(x))})
  dataExpNorm           = normalize(dataExp)
  #plot Exp density
  if (plotflag == 1){
    outPdfName = paste(outputcoeff, genename,".pdf",sep="")
    if(file.exists(outPdfName)){
      pdf(paste(outputcoeff, genename,"_cnv.pdf",sep=""))
    }else{
      pdf(outPdfName)
    }
    plot(density(dataExp),main = "Gene Expression Density")
    plot(density(dataExpNorm), main = "Gene Expression Density after normalization")
  }

  #---------loading snp data
  dataSnp             = read.table(inputsnp,header =T)
  rownames(dataSnp)   = dataSnp[,2]
  dataSnp             = dataSnp[,-c(1:2)]
  #snp data QC:1. delete all 0 snp
  dataSnp             = dataSnp[rowSums(dataSnp)>0,]

  if (nrow(dataSnp) > 0 ){
      ##load cnv data
      dataCnv             = read.table(inputcnv,header =T)
      rownames(dataCnv)   = dataCnv[,2]
      dataCnv             = dataCnv[,-c(1:2)]

      ##------
      data_merge  = t(rbind(dataExpNorm,dataCnv,dataSnp))
    }
  return(list(dataSnp=dataSnp,
              dataExp=dataExp,
              dataCnv=dataCnv,
              data_merge=data_merge,
              genename = genename))
}
print("loading data")   
allData   = getData(inputexp,inputsnp,inputcnv,plotflag)
dataSnp   = allData$dataSnp
dataCnv   = allData$dataCnv
dataExp   = allData$dataExp
data_merge= allData$data_merge
genename  = as.character(allData$genename)
outputFtest = paste(outputFtest,genename,".txt",sep="")
print(paste("output ftest file:",outputFtest ))
#########-----------model------------

if(nrow(dataSnp > 0)){

    ##---------fitting cnv + snp
    grpRegPermuFile     = paste("grplasso_",genename,".RData",sep="")
    if( file.exists(grpRegPermuFile) ){
      print("loading permutation RData")
      load(grpRegPermuFile)
    }else{   
      ##-----------------Grouping_of_variables
      group = groupVars(dataSnp,dataCnv,plotflag)
      require(grpreg)    
      print("Doing regression...")
      nperm     = 1000
      fitpermut = regfitPermu(data_merge, group, nperm, plotflag)
      ##saving data
      save(fitpermut,group,data_merge,
           file=paste("grplasso_",genename,".RData",sep=""))
      ##output
      print("writing output...")
      outputcoeff = paste(outputcoeff,genename,".txt",sep="")
      write.table(t(as.matrix(c(paste("gene","RSS","npermu","pvalue",sep=":"),paste(genename,fitpermut$RSS, nperm, fitpermut$pvalue,sep=":")))),
                  outputcoeff,
                  quote=F,col.names=F,sep="\t",row.names = F)
      write.table(as.matrix(sort(fitpermut$beta)),
                  outputcoeff,
                  append=T,
                  col.names=F,quote=F,sep="\t")
    }
   ##---------fitting cnv
   if (nrow(dataCnv) ==1){
      print("fitting cnv...")    
      X = data_merge[,2]
      y = sigmoid(data_merge[,1]) 
      cnvFit = fitCnv(X,y,plotflag)
      
      ##F-test of residual from cnv and residual from grplasso
      print("doing F test...")
      p1 = nrow(dataCnv)
      p2 = nrow(dataSnp) + p1
      nPoints = ncol(dataSnp)
      regFtest = myFtest(cnvFit$residuals,fitpermut$fit_residual,p1,p2,nPoints)
      RSScnv = sum(cnvFit$residuals^2)
      RSSgrp = fitpermut$RSS
      print("output Ftest result")
      write.table(t(as.matrix(c(paste("gene","RSScnv", "RSSgrp", "fStatistic", "pvalue",sep=":"),
                                paste(genename, RSScnv, RSSgrp, regFtest$f, regFtest$pval,sep=":")))),
                  outputFtest,
                  quote=F,col.names=F,sep="\t",row.names = F)
    }else{
      print("No cnv around...")
      residual_cnv = rep(0,length(dataExp))
    }  
    
}else {
  write.table(t(as.matrix(c(paste("gene","RSS","npermu","pvalue",sep=":"),paste(genename,'NA', 'NA', 'NA',sep=":")))),
                  outputcoeff,
                  quote=F,col.names=F,sep="\t",row.names = F)
}