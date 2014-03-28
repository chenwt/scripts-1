#!/usr/bin/Rscript
#Author: Jing He
#Date:24 Oct,2013 
#Last Updated:
#COMMENTS: need edgeR installed; 
#input: <string:path you wnat your results to be> 
# 		  <string:name of your design file(4 cols, tab delimite:example)
#		    <string:name of count matrix file>
# 		  <string:name of your output files>
#output: <file:pdf of 2 plots> <file: txt of differetial expresssed genes>
####TODO: need more development

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  rootd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
}

args = getArgs()
usage = "Usage: Rscript step1-2_expVoomNormed.r --tumor <tumor.mat file> --normal <normal.mat file> "
example = "Example: "
# if(length(args) < 3 || is.null(args)){
#   print(usage)
#   print(example)
#   print(args)
#   stop("Input parameter error!")
# }else{
#   print(args)
# }

CDT = as.character(paste(unlist(strsplit(date(),split=" "))[c(2,3,5)],collapse="-"))
setwd(system("pwd",intern=T))
cwd         = getwd()
# tumor  = args['tumor'] 
# normal = args['normal']
# outExp = paste(tumor, "_", CDT, ".voomNormed.matrix", sep="")
print(paste("current working directory:",cwd))


##-----test
# cwd = getwd()
tumor  =  "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix"
normal =  "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_normal_Mar-21-2014.matrix"
outExp = paste(tumor, "_", CDT, ".voomNormed.matrix", sep="")
outExpN = paste(normal, "_", CDT, ".voomNormed.matrix", sep="")

##----------------------------
getData = function(file,type="T"){
  data = read.delim(file,header=T)
  gene = unique(unlist(data[,1]))
  data = data[-16274,]
  data = data[,-1]
  rownames(data)=gene
#   colnames(data) = vapply(colnames(data),FUN=function(x){substr(x,start=9,stop=16)},'a')
  sample = colnames(data)
  design = rep(type,ncol(data))
  names(design) = sample
  return(list(data=data,design=design,gene=gene))
}

##--------#load data

dataT = getData(tumor,type='tumor')
dataN = getData(normal,type='normal')
cntSampleT = ncol(dataT$data)
cntSampleN = ncol(dataN$data)
cntGene = nrow(dataT$data)
dataMat = cbind(dataT$data,dataN$data)
gene = dataT$gene
row.names(dataMat) = gene

design = c(dataT$design,dataN$design)
condition <- factor( design )

##-------------voom transformation
require(limma)
designMat = model.matrix(~design)
dataMatVoom = voom(as.matrix(dataMat),designMat,plot=TRUE)

dataMatVoomTumor = dataMatVoom$E[,which(as.character(design)=="tumor")]
dataMatVoomNormal = dataMatVoom$E[,which(as.character(design)=="normal")]

##---------output
dev.off()
write.table(round(dataMatVoomTumor,digits=4),
            outExp,sep="\t",col.names=T,row.names=T,quote=F)

write.table(round(dataMatVoomNormal,digits=4),
            outExpN,sep="\t",col.names=T,row.names=T,quote=F)

#--------------------------------------------
