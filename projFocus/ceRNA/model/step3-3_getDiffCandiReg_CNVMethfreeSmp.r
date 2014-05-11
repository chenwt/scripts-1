#!/usr/bin/Rscript
#Author: Jing He
#COMMENTS: after get CNV, methylation free sample  list, for each candidate dirver, 
# get thos candidate drivers which differential expressed when compare mutated to non-mutated samples
#input: <string:path you wnat your results to be> 
# 		  <string:name of your design file(4 cols, tab delimite:example)
#		    <string:name of count matrix file>
#output: <file: txt of differetial expresssed genes>


sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/myR/jingGraphic.R")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  rootd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/myR/jingGraphic.R")
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
# 
# setwd(system("pwd",intern=T))
# cwd         = getwd()
# tumor  = args['tumor'] 
# normal = args['normal']
# output = args['out'] 
# gslist = args['gslist']
# print(paste("current working directory:",cwd))
# 

##-----test
cwd = getwd()
tumor  =  "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
output =  jxy("brca_mutCandiReg_diff",substr(Sys.time(),1,10))
gslist = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step3_funcMutKegReg/regIntSmp/gslist_Apr-23-2014_CnvMethFree.nonzeor"
##-----test

CDT = gsub("-","",substr(Sys.time(),1,10))
output = paste(gslist,".deg_", CDT,".txt",sep="")
# outExp = paste(tumor,"_Voomed_DEGexp_", CDT,".rda",sep="")
# outputPDF = paste(reportDir,"/",output,".pdf",sep="")
##----------------------------
getData = function(file,type="T",glist){
  gene = unlist(read.table(glist)[,1])
  data = read.delim(file,header=T)
  gene = gene[ gene %in% data[,1] ]
  data = data[data[,1] %in% gene,-1]
  rownames(data)=gene
  sample = colnames(data)
  design = rep(type,ncol(data))
  names(design) = sample
  return(list(data=data,design=design,gene=gene))
}

##--------#load data
dataT = read.table(tumor,sep="\t",header=T)
cntSampleT = ncol(dataT$data)
cntGene = nrow(dataT$data)
gene = dataT
expM = dataT
allSmpExp = colnames(expM)
allexpGene = rownames(expM)
topleft(expM)

mutfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix"
mutM = read.table(mutfile, header = F)
mutSmp = vapply(mutM[1,][-1],as.character,'a')
if (nchar(mutSmp[1]) > 9)  mutSmp = vapply(mutSmp, FUN = function(x){gsub("-",".",substr(x,6,16))},'a')
colnames(mutM) = mutSmp
rownames(mutM) = mutM[,1]
nlast = ncol(mutM)
mutM = mutM[,-c(1,nlast)][-1,]
allmutGene = rownames(mutM)
# require(useful)
# topright(mutM)
numMut = nrow(mutM)

##load intact sample information
con  <- file(gslist, open = "r")
intactSmpL = list()
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  templine = unlist(strsplit(oneLine,"\t"))
  intactSmpL[templine[1]] = strsplit(templine[2],";")
} 
close(con)

calZscore = function(x, y){
#   x = x[1,]
#   y = y[1,]
  return((sum(x)/length(x) - sum(y)/length(y)) / sd(y))
}

out = as.data.frame(matrix(NA,nrow=numMut,ncol=7))
colnames(out) = c('regGeneName', "numIntSmp", "numMutIntSmp","numMut", "numNonMut", "zscore","pvalue")
for (i in 1:numMut ){
#   i = 4
  mutV = vapply(mutM[i,],function(x){as.numeric(as.character(x))},0)
  geneCrt = allmutGene[i]
  tmpIntSmp = intactSmpL[[geneCrt]]
  tmpExp  = expM[geneCrt,intersect(tmpIntSmp,allSmpExp)]
  tmpMutSmp = intersect(intersect(names(mutV[which(mutV!=0)]),tmpIntSmp),allSmpExp)
  if (length(grep(paste("^",geneCrt,"$",sep=""),allexpGene)) > 0 && length(tmpMutSmp) > 0 ){
#       tmpExp = vapply(tmpExp,function(x){as.numeric(as.character(x))}, 1.0)
      mutIdx  = vapply(tmpMutSmp,FUN=function(x){grep(x,names(tmpExp))},1)
      expMutTempV = as.vector(tmpExp[mutIdx])
      expNonMutTempV = as.vector(tmpExp[-mutIdx])
      tmpZs = calZscore(expMutTempV, expNonMutTempV)
      out[i,] = c(geneCrt,length(tmpIntSmp) , length(tmpMutSmp), 
                  length(expMutTempV), length(expNonMutTempV),
                  tmpZs, 2*pnorm(-abs(tmpZs)) )
      next
  }else{
    out[i,] = c(geneCrt,length(tmpIntSmp), length(tmpMutSmp),NA,NA, NA, NA)
    print(paste(geneCrt,"gene not in exp data"))
    next
  }
}

length(na.omit(out[order(as.numeric(out$pvalue)),5]))
out1 = out[order(as.numeric(out$pvalue)),]
write.table(out1, file = output, quote = F, sep = "\t", row.names=F)
out2 = na.omit(out1)
out2$pValAdj = p.adjust(out2$pvalue,method='fdr')
out2  = out2[order(as.numeric(out2$pValAdj)),]
write.table(out2, file = jxy(output,".adj"), quote = F, sep = "\t", row.names=F)
