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
  setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression/")
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  figd = "/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/report/figure/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/expression")
  rootd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  figd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}

# args = getArgs()
# usage = "Usage: Rscript step1-2_DEG_UCceRNET_edgeR.r --tumor <tumor.mat file> --normal <normal.mat file>  --genelist <geneSample.list>  "
# example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterGrplasso.r --file grplasso_coeff --cut 0.05 --out gcGenes_GeneVarNet"
# if(length(args) < 3 || is.null(args)){
#   print(usage)
#   print(example)
#   print(args)
#   stop("Input parameter error!")
# }else{
#   print(args)
# }

setwd(system("pwd",intern=T))
# cwd         = getwd()
# tumor  = args['tumor'] 
# normal = args['normal']
# output = args['out'] 
# print(paste("current working directory:",cwd))


##-----test
cwd = getwd()
tumor  =  "brca_exp_level3_02042014.mat"
normal =  "brca_expNormal_level3_02042014.mat"
# output =  "brca_ucCeRNETCancerGeneDEG_edgeR_02042014"
crtDate = substr(Sys.time(),1,10)
outExp =  paste(tumor,"_voomed_",crtDate,".mat",sep="")
reportDir = paste(rootd, "/report/topDown_02042014/fig",sep="")
outputPDF = paste(reportDir,"/",outExp,".pdf",sep="")
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
# dataMatVoomNormal = dataMatVoom$E[,which(as.character(design)=="normal")]

##---------output
# idxOutGene = samplePerGene$gene %in% setdiff(tarGene,failGene)
# out  =samplePerGene[idxOutGene,]
# rownames(out) = samplePerGene$gene[idxOutGene]
# write.table(samplePerGene[idxOutGene,],outGeneSample ,sep="\t",col.names=T,row.names=F,quote=F)
# 
# expTumor = dataMatVoomTumor[rownames(dataMatVoomTumor) %in% out$gene,]
# expNormal =dataMatVoomNormal[rownames(dataMatVoomNormal) %in% out$gene,]
write.table(dataMatVoomTumor,outExp,sep="\t",col.names=T,row.names=T,quote=F)
#--------------------------------------------
