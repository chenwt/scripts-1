##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <string: regulator gene name > 
#output: <file: > 
#Usage: Rscript grplassoSNP.r input
#Description:  do lasso to each cancer gene; call grplassoSNP.r for regression of each of it's target  
##           plot expression for target and all it's regulators
example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript 
  /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-2_targetLasso.r --gene ESR1 --out ESR1_grplasso "
usage = "Usage: Rscript step2-2_targetLasso.r  --gene <gene name> "
ERR = "ERROR here: "

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  rootd= "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  source(paste(rootd, "/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
  source(paste(rootd, "/scripts/projFocus/ceRNA/model/grpLassoFunctions.R",sep=""))
  source(paste(rootd, "/scripts/myR/jingGraphic.R",sep=""))
  setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/")

}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/grpLassoFunctions.R")
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/myR/jingGraphic.R")
  print("working from Linux")
  rootd= "/ifs/home/c2b2/ac_lab/jh3283/"
  setwd("/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/")
}

#----getting command line parameters
# args = getArgs()
# if(length(args) < 1 || is.null(args)){
#   print(paste(usage,example,sep="\n"))
#   print(args)
#   stop(paste(ERR,"wrong input parameter!"))
# }

setwd(system("pwd",intern=T))
# cwd         = getwd()
# gene        = args['gene']
# out         = args['out']
# output      = paste(cwd,"/",args['out'], sep="")

#--test
cwd         = getwd()
gene        = "ESR1"
out         = paste(gene,"grplasso",sep="")
output      = paste(cwd,"/",out, sep="") 

##-------------prepareData
# cmd = paste(rootd,"/scripts/projFocus/ceRNA/model/getGeneData.sh  ", gene, sep="")
# system(cmd,intern=T, wait=T)
setwd(paste(cwd,"/temp-",gene,sep=""))

##get samples
samples = vapply(unlist(strsplit(unique(read.table("sample",stringsAsFactors=F)[,3]),";")),subStr1To19,'a')
###---target gene
expTarget = formatData(inputFile=paste(gene,"_exp",sep=""),t="exp")
expTarget = subset(expTarget,select=intersect(samples,colnames(expTarget)))
expTarget = t(as.matrix(expTarget[,order(expTarget[1,])]))
# expTarget = t(as.matrix(expTarget[,colnames(expTarget) %in% samples]))
rownames(expTarget) = gene

##get all regulators
regulators = read.table("stat.txt",stringsAsFactors=FALSE,skip=1,header=T)
regulators = regulators[regulators[,4]>0,]
regulators = unique(vapply(regulators[,1],as.character,'a'))
numReg  = length(regulators)

##----get expression of all regulators
expRegulators = as.data.frame(matrix(0,nrow=numReg,ncol=length(expTarget)))
so = colnames(expTarget)
for (i in 1:numReg){
  tempReg = regulators[i]
  tempExp  = formatData(paste(tempReg,"_exp",sep=""),t="exp")
#   xx = intersect(samples,colnames(tempExp))
  expRegulators[i,] = subset(tempExp,select=so) 
}


##---diagnostic plot
#  par(mar=c(4,4,1,0))
#  hist(regulators[,3],breaks=14,col="blue",
#      xlab=colnames(regulators)[3],ylab="Number of SNP",
#      main="")

###-------------plot
drawVectorHeat = function(z, pos = 0, size = 0.5) {
  myCol = colorRampPalette(c("blue","white","red"))(256)
  z = val2col(z,col=myCol)
  for (i in 1:length(z)){
    rect(i-1,ybottom=pos ,xright=i,ytop=(pos+ size),col=z[i],border=z[i])  
  }
}
addPIDLabel = function(z){
  colLabels = colnames(z)
  axis(1,at=seq(0.5,length(z),1),labels=colLabels,cex.axis=0.4,
       las=2, lwd=0.3, tck = -0.01, col.ticks="gray")
}
addGeneLable = function(label,pos){
  axis(2,at=pos,cex.axis=0.5, 
       labels = label,las=2,lwd=0.3, tck = -0.01,  col.ticks="gray")
}

system("mkdir fig")
pdf("fig/expression.pdf")
par(mar=c(4,6,1,0))
blankPlot(maxX=length(expTarget),maxY=length(regulators))
drawVectorHeat(expTarget,pos=0,size=2)
addGeneLable(pos=1,label=gene)
addPIDLabel(expTarget)
for (i in 1:numReg){
  size = 1
  drawVectorHeat(unlist(expRegulators[i,]),pos=(i * size) + 2,size=size)
  addGeneLable(pos= i*size - size/2 + 3,label=regulators[i])
}
dev.off()
# 
# require(gplots)
# heatmap.2(as.matrix(expRegulators),col=colorRampPalette(c("blue","white","red"))(256),trace = "none")
# expCol = colorRampPalette(c("blue", "red"))(max(regulators[,1]))
# cnvCol = colorRampPalette(c("yellow","red"))(max(regulators[,2]) + 1)
# snpCol = colorRampPalette(c("white","black")[max(regulators[,3]) + 1])
# somCol = "white"
# dataColor = as.data.frame(matrix("",nrow=numR,ncol=numC))
# for (i in 1:numR){
#   dataColor[i,1] = expCol[data[i,1]]
#   dataColor[i,2] = cnvCol[exp(data[i,2]) + 1]
#   dataColor[i,3] = snpCol[data[i,3] + 1]
#   dataColor[i,4] = somCol
# }
# 
# 
# ##----group all regulators into different analysis type based on mutational spectrum
# # 
# # colnames(regulators) = regulators[1,]
# # numGene = nrow(regulators)
# # drawTable(regulators)
# # for ( i in 1:numGene){
# #   i=154
# #   temp = regulators[i,]
# #   tempGene = as.character(unlist(temp[1]))
# #   ##focus on somatic mutation
# #   if (temp$CNV > 0  && temp$SNP > 0  && temp$SOM > 0) 
# #       type = 1
# #       cmd = paste("grep -w ", tempGene, " ") #<<==== HERE
# # }
# 
