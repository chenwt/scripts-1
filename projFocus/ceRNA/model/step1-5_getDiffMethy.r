##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file1: tumor methlation matrix> <file2: normal methylation matrix> 
#output: <file: tumor sample with methylation relative to population mean>
#Description: this file was created for projFocus, ceRNA, used in step1 to eliminate tumor samples
#TODO: 

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth/")
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  figd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth")
  rootd = "/ifs/scratch/c2b2/ac_lab/jh3283/"
  figd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}

args = getArgs()
usage = "Usage: Rscript step1-5_getDiffMethy.r --tumor <tumor mat> --normal <normal mat> --out <output figure to fig>  "
example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript  step1-5_getDiffMethy.r brca_methNormal_level3_02072014.mat brca_methTumor_level3_02072014.mat"
if(length(args) < 2 || is.null(args)){
  print(usage)
  print(example)
  print(args)
  stop("Input parameter error!")
}else{
  print(args)
}

###-----functions------

fillNA = function(var){
  var[is.na(var)] <- mean(var, na.rm = TRUE)
  return(var)
}

betaToM = function(beta){
  max = 50
  min = -50
  b = as.numeric(beta)
  m = ifelse(b >= 1,10,ifelse(beta <= 0,-10,log2(beta /(1 - beta)) ) )
  return(m)
}

#-----end--func
setwd(system("pwd",intern=T))
cwd         = getwd()
normal      = args['normal']
tumor       = args['tumor']
# outfile     = args['out']

require('sva')
library(charm)


dataT = read.delim(tumor,header=T)
row.names(dataT) = dataT$barcode
dataT = apply(dataT[,-1],c(1,2),as.numeric)
dataT = dataT[,-ncol(dataT)]

dataN = read.delim(normal,header=T)
row.names(dataN) = dataN$barcode
dataN = apply(dataN[,-1],c(1,2),as.numeric)
dataN = dataN[,-ncol(dataN)]

##remove outlier normal samples
meansN = colSums(dataN)/nrow(dataN)
dataN = dataN[,which((meansN - mean(meansN)) / sd(meansN) < 2)]
dataAll = cbind(dataT,dataN)
dataAll = apply(dataAll,2,fillNA)

##----normalization
batch = as.matrix(c(rep("methylation27k",46),rep("methylation450k",(ncol(dataAll)-46))))
mod = model.matrix(~as.factor(c(rep("tumor",ncol(dataT)),rep("normal",ncol(dataN)) ) ))
indexNA27k = which(rowSums(dataT[,1:46])==0,arr.ind=T)
dataAll.norm = ComBat(dataAll[-indexNA27k,],batch=batch,mod=mod)

dataMAll.norm = apply(dataAll.norm,c(1,2),betaToM)
date = Sys.Date()
# pdf(paste(figd,outfile,Sys.Date(),".pdf",sep=""))
# plotValue(dataAll[,1:117],dataAll[,118:205])
# plotValue(dataMAll.norm[,1:117],dataMAll.norm[,118:205])
# dev.off()
# 

###-------------statistical testing
dataT.norm = dataMAll.norm[,1:117]
dataN.norm = dataMAll.norm[,118:205]
genMethMat = function(dataT.norm, dataN.norm, zcut = 2.5){  
  meansNormal = rowSums(dataN.norm) /ncol(dataN.norm)
  dataT.res = matrix(0,ncol=ncol(dataT.norm),nrow=nrow(dataT.norm))
  colnames(dataT.res) = colnames(dataT.norm)
  rownames(dataT.res) = rownames(dataT.norm)
  for (n in 1:ncol(dataT.norm)){
    dataT.res[,n] =  ifelse((dataT.norm[,n] - meansNormal)/sd(meansNormal) > zcut,1,0)
  }
  return(dataT.res)
}

result = genMethMat(dataT.norm,dataN.norm,zcut = 2.5)
out = paste(tumor,"diffMeth.mat",sep="_")
write.table(result,out,sep="\t",quote=F)


