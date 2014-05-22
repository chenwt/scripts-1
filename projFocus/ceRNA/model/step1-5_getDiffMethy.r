##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file1: tumor methlation matrix> <file2: normal methylation matrix> 
#output: <file: tumor sample with methylation relative to population mean>
#Description: this file was created for projFocus, ceRNA, used in step1 to eliminate tumor samples
#TODO: 

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth/")
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  figd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014//fig/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth")
  rootd = "/ifs/scratch/c2b2/ac_lab/jh3283/"
  figd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}

# args = getArgs()
# usage = "Usage: Rscript step1-5_getDiffMethy.r --tumor <meth_tumor.matrix> --normal <meth_normal.matrix> --out <meth_diff.matrix> 
# --t27 315 --n27 27 "
# example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript 
#     --tumor 
#     --normal
#     "
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
CDT = paste(unlist(strsplit(date(),split=" "))[c(2,3,5)],collapse="-")
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

###---------funcEnd-------
setwd(system("pwd",intern=T))
# cwd         = getwd()
# tumor       = args['tumor']
# normal      = args['tumor']
# numTum27    = as.numeric(args['t27'])
# numNorm27   = as.numeric(args['n27'])
print(paste("current working directory:",cwd))

## source("http://bioconductor.org/biocLite.R")
## biocLite("charm")
## source("http://bioconductor.org/biocLite.R")
## biocLite("sva")
# library(charm)
##-----test---
normal = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth/brca_meth_l3_normal_Mar-23-2014.matrix"
tumor =  "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/meth/brca_meth_l3_tumor_Mar-24-2014.matrix"
numTum27 = 315
numNorm27 = 27
##----testend----


dataT = read.delim(tumor,header=T)
row.names(dataT) = paste(dataT[,1], dataT[,2],sep=":")
dataT = apply(dataT[,-c(1,2)],c(1,2),as.numeric)
dataT = dataT[,-ncol(dataT)]
dataN = read.delim(normal,header=T)
row.names(dataN) = paste(dataN[,1], dataN[,2],sep=":")
dataN = apply(dataN[,-c(1,2)],c(1,2),as.numeric)
dataN = dataN[,-ncol(dataN)]

# pdf(jxy(figd, "/scatterplot_meth_", CDT, ".pdf"))
# plotValue(dataT,dataN)
# dev.off()
##remove outlier normal samples
meansN = colSums(dataN)/nrow(dataN)
dataN = dataN[,which((meansN - mean(meansN)) / sd(meansN) < 2)]

numT = ncol(dataT)
numN = ncol(dataN)

dataAll = cbind(dataT,dataN)
dataAll = apply(dataAll,2,fillNA)

##----normalization
batch = as.matrix(c(rep("methylation27k",numTum27),rep("methylation450k",(numT - numTum27)),
                    rep("methylation27k",numNorm27),rep("methylation450k",(numN - numNorm27)) ))

mod = model.matrix(~as.factor(c(rep("tumor",numT),rep("normal",numN ) ) ))
# indexNA27k = which(rowSums(dataT[,c(1:numTum27, numT+1:numT+numNorm27])==0,arr.ind=T)
require('sva')
dataAll.norm = ComBat(dataAll, batch = batch, mod = mod)

dataMAll.norm = apply(dataAll.norm,c(1,2),betaToM)

###-------------statistical testing
dataT.norm = dataMAll.norm[,1:numT]
dataN.norm = dataMAll.norm[,(numT + 1):(numT + numN) ]

result = genMethMat(dataT.norm,dataN.norm,zcut = 2.0)
outGenes = unique((unlist(sapply(rownames(result),FUN=function(x){unlist(strsplit(
                          unlist(strsplit(x,":"))[2], ";"))}))))

probeGenes = rownames(dataAll.norm)
outResult = as.data.frame(matrix(0,nrow=length(outGenes), ncol = ncol(result)))
colnames(outResult) = colnames(result)
rownames(outResult) = outGenes
for (i in 1: length(outGenes)){
  ptn = jxy(":", outGenes[i], "$", "|:", outGenes[i], "|;", outGenes[i], "|", outGenes[i],"$")
  tempIdx = grep(ptn,probeGenes)
  if (length(tempIdx) > 1){
    outResult[i,] = apply(result[tempIdx,],2,function(x){x[which.max(abs(x))]})
  }else if (length(tempIdx) == 1){
    outResult[i,] = result[tempIdx,]
  }
}
out = paste(tumor, "diffMeth.matrix",sep="_")
write.table(outResult,out,sep="\t",quote=F)

# pdf(jxy(figd,"/scatterplot_meth_level3_normMvalue_",CDT, ".pdf"))
# plotValue(dataAll.norm[,1:numT],dataAll.norm[,(numT + 1):(numT + numN)])
# plotValue(dataMAll.norm[,1:numT],dataMAll.norm[,(numT + 1):(numT + numN)])
# barplot(rowSums(result),width=4,col="blue")
# hist(table(rowSums(result)))
# image(result)
# dev.off()

###------------------------ descriptive examing data
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
# plotValue(dataT,dataN)
# plotValue(dataMT,dataMN)
plotDensity = function(dataT,dataN, n =5){
  layout(matrix(1:4,nrow=2,byrow=T))
  for (i in n:(n+3)){
    print(plot(density(dataT[i,]),main = rownames(dataT)[i],col="red"))
    print(lines(density(dataN[i,]),main = rownames(dataN)[i],col="blue"))
  }
}
# plotDensity(dataT,dataN)
# require(limma)
# dataAll.normMedian = normalizeMedianValues(dataAll)
# boxplot(dataAll.norm)
# heatmap.2(na.omit(dataAll.norm), trace= "none", main="All meth", col=mycol)
# dataT.norm = normQ(dataT)
# dataN.norm = normQ(dataN)
# plotValue(dataT.norm,dataN.norm)
# plotValue(dataMT,dataMN)
# plotDensity(dataAll.norm[,1:96],dataAll.norm[,97:213])
# plotDensity(dataAll.normMedian[,1:96],dataAll.normMedian[,97:213])
# dataAll.normBtwArray = normalizeBetweenArrays(dataAll)
# plotDensity(dataAll.normBtwArray[,1:96],dataAll.normBtwArray[,97:213])
# plotValue(dataMT.scale,dataMN.scale)
# 
