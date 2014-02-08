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
# 
# args = getArgs()
# usage = "Usage: Rscript bridegCeRAN.r --file <gene.list>  "
# example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterGrplasso.r --file grplasso_coeff --cut 0.05 --out gcGenes_GeneVarNet"
# if(length(args) < 3 || is.null(args)){
#   print(usage)
#   print(example)
#   print(args)
#   stop("Input parameter error!")
# }else{
#   print(args)
# }
# 

###-----functions------


fillNA = function(var){
  var[is.na(var)] <- mean(var, na.rm = TRUE)
  return(var)
}
setwd(system("pwd",intern=T))
cwd         = getwd()
# filepref    = args['file']
# cutoff      = as.numeric(args['cut'])
# output      = paste(args['out'],".rda",sep="")
# print(paste("current working directory:",cwd))
# source("http://bioconductor.org/biocLite.R")
# biocLite("charm")
source("http://bioconductor.org/biocLite.R")
biocLite("sva")
require('sva')
library(charm)
normal = "brca_methNormal_level3_02072014.mat"
tumor =  "brca_methTumor_level3_02072014.mat"

dataT = read.delim(tumor,header=T)
row.names(dataT) = dataT$barcode
dataT = apply(dataT[,-1],c(1,2),as.numeric)
dataT = dataT[,-ncol(dataT)]
dataN = read.delim(normal,header=T)
row.names(dataN) = dataN$barcode
dataN = apply(dataN[,-1],c(1,2),as.numeric)
dataN = dataN[,-ncol(dataN)]

# dataT = apply(dataT, c(1,2),betaToM)
# dataN = apply(dataN, c(1,2),betaToM)
plotValue(dataT,dataN)
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



betaToM = function(beta){
  max = 50
  min = -50
  b = as.numeric(beta)
  m = ifelse(b >= 1,10,ifelse(beta <= 0,-10,log2(beta /(1 - beta)) ) )
  return(m)
}
dataMAll.norm = apply(dataAll.norm,c(1,2),betaToM)

pdf(paste(figd,"meth_level3_norm_02072014.pdf",sep=""))
plotValue(dataAll[,1:117],dataAll[,118:205])
plotValue(dataMAll.norm[,1:117],dataMAll.norm[,118:205])
dev.off()


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
# barplot(rowSums(result),width=4,col="blue")
# hist(table(rowSums(result)))
# image(result)
out = paste(tumor,"diffMeth.mat",sep="_")
write.table(result,out,sep="\t",quote=F)

### descriptive examing data
plotValue = function(dataT,dataN){
    layout(matrix(1:4,nrow=2,byrow=T))
    boxplot(dataT[,-1], main= "tumor meth" , xlab= "sample")
    boxplot(dataN[,-1], main= " normal meth",xlab= "sample")
    cut = 47
    boxplot(dataT[,1:47],main="tumor meth27",xlab= "sample")
    boxplot(dataT[,48:ncol(dataT)],main="tumor meth450",xlab= "sample")
}
plotGene = function(dataT,dataN){
  require('gplots')
  mycol = bluered(256)
  heatmap.2(na.omit(dataT), trace= "none", main="tumor meth", col=mycol)
  heatmap.2(na.omit(dataT[,1:47]), trace= "none", main="tumor meth27", col=mycol)
  heatmap.2(na.omit(dataT[,48:96]), trace= "none", main="tumor meth450", col=mycol)
  heatmap.2(na.omit(dataN), trace= "none", main="normal meth", col=mycol)
}
plotValue(dataT,dataN)
plotValue(dataMT,dataMN)
plotDensity = function(dataT,dataN, n =5){
  layout(matrix(1:4,nrow=2,byrow=T))
  for (i in n:(n+3)){
    print(plot(density(dataT[i,]),main = rownames(dataT)[i],col="red"))
    print(lines(density(dataN[i,]),main = rownames(dataN)[i],col="blue"))
  }
}
plotDensity(dataT,dataN)
require(limma)
dataAll.normMedian = normalizeMedianValues(dataAll)
boxplot(dataAll.norm)
heatmap.2(na.omit(dataAll.norm), trace= "none", main="All meth", col=mycol)
dataT.norm = normQ(dataT)
dataN.norm = normQ(dataN)
plotValue(dataT.norm,dataN.norm)
plotValue(dataMT,dataMN)
plotDensity(dataAll.norm[,1:96],dataAll.norm[,97:213])
plotDensity(dataAll.normMedian[,1:96],dataAll.normMedian[,97:213])
dataAll.normBtwArray = normalizeBetweenArrays(dataAll)
plotDensity(dataAll.normBtwArray[,1:96],dataAll.normBtwArray[,97:213])
plotValue(dataMT.scale,dataMN.scale)

