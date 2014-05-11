#!/usr/bin/Rscript
#Author: Jing He
#COMMENTS: after get CNV, methylation free sample  list, for each candidate dirver, 
# get thos candidate drivers which differential expressed when compare mutated to non-mutated samples
#input:1. group lasso result keyRegsFile
#        2. geneSample list of intact samples for each gene
#        3. voom normalized expression file for all targe genes

#output: enrichment score for  ceRNA regulator mutaion with its target genes


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
CDT = gsub("-","",substr(Sys.time(),1,10))

setwd("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step3_funcMutKegReg/tarRegMutAssocation/")
cwd = getwd()
gslistf = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethSomFree.10smapMore.deg_20140325.txt.10more.hasReg.list"
tarDriverf = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/summary/target_ceRNADriver_0.01"
expTufilef = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
mutfilef = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/sigMut/step2_mutKeyReg/kegRegs_Apr-18-2014.summary.driverRegs.list.uniq.mut.matrix"
figd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Apr2014/fig/"
output = paste(cwd,"/","diff_TarRegMut_tarIntSmp_",CDT,".tsv",sep="")

##-----test

##--------#load data

con  <- file(gslistf, open = "r")
intactSmpL = list()
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  templine = unlist(strsplit(oneLine,"\t"))
  intactSmpL[templine[1]] = strsplit(templine[2],";")
} 
close(con)

con  <- file(tarDriverf, open = "r")
tDrCeRNAL = list()
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  templine = unlist(strsplit(oneLine,"\t"))
  tDrCeRNAL[templine[1]] = strsplit(templine[2],";")
} 
close(con)


expD = read.table(expTufilef, header=T, sep="\t")
tarExpSmpV = colnames(expD)

require(useful)


mutM = read.table(mutfilef, header = F)
mutSmp = vapply(mutM[1,][-1],as.character,'a')
if (nchar(mutSmp[1]) > 9)  mutSmp = vapply(mutSmp, FUN = function(x){gsub("-",".",substr(x,6,16))},'a')
allmutGene = mutM[,1]
mutM = mutM[,-1][-1,]
colnames(mutM) = mutSmp
rownames(mutM) = allmutGene[-1]
nlast = ncol(mutM)



numTarget = length(tDrCeRNAL)
allTargets = names(tDrCeRNAL)

calZscore = function(x, y){
  return((sum(x)/length(x) - sum(y)/length(y)) / sd(y))
}

myZtest = function(m, e){
  mm = vapply(m,FUN=function(x){as.numeric(as.character(x))},1)
  smpm = names(m)
  smpmm1 = smpm[which(mm!=0)]
  if (length(smpmm1) > 0) {
    smpmm2 = setdiff(smpm,smpmm1) 
    tmpZs = calZscore(unlist(e[smpmm1]), unlist(e[smpmm2]))
    tmpPv = 2*pnorm(-abs(tmpZs))
    out = c(target = tgene,ceRNA = row.names(m), numMut= length(smpmm1), numNonMut = length(smpmm2),zscore = tmpZs,pval = tmpPv)
    
    if (tmpPv < 0.1 | is.na(tmpPv)) {
      print(paste(tgene,row.names(m),tmpZs,tmpPv))
    }
    return(out)
    
  }else{
    return(rep(NA,6))
  }
}


###-----for each 
outD = as.data.frame(matrix(NA,ncol=6))
for(j in 1:length(allTargets)){
# for(j in 1:10){
  tgene = allTargets[j]  
  mutIntExpSmp = intersect(intactSmpL[[tgene]],intersect(tarExpSmpV,mutSmp))
  if (length(mutIntExpSmp) > 0) {
    tExpD = subset(expD,select=mutIntExpSmp)[tgene,]
    regMut = subset(mutM,select=mutIntExpSmp)[intersect(allmutGene,tDrCeRNAL[[tgene]]),]
    if (nrow(regMut) > 0 & nrow(na.omit(tExpD))){
      for (i in 1:nrow(regMut)){
       crtOut = myZtest(regMut[i,], tExpD)    
       outD = rbind(outD,crtOut)
      }
    }
  }
}
outD = na.omit(outD)
colnames(outD) = c("target","ceRNA","numMut", "numNonMut","zscore","pval")
outD$pvalAdj = p.adjust(outD$pval,method='fdr')


out1 = outD[order(as.numeric(outD$pval)),]
head(out1)
write.table(out1, file = output, quote = F, sep = "\t", row.names=F)


###---agreagate---
outD = as.data.frame(matrix(NA,ncol=6))
for(j in 1:length(allTargets)){
#   for(j in 1:10){
#     j = 10
  tgene = allTargets[j]  
  mutIntExpSmp = intersect(intactSmpL[[tgene]],intersect(tarExpSmpV,mutSmp))
  if (length(mutIntExpSmp) > 0) {
    tExpD = subset(expD,select=mutIntExpSmp)[tgene,]
    regMut = subset(mutM,select=mutIntExpSmp)[intersect(allmutGene,tDrCeRNAL[[tgene]]),]
    if( nrow(regMut) > 0 & nrow(na.omit(tExpD)) > 0 ){
      regMut = apply(regMut,c(1,2),function(x){as.numeric(as.character(x))})
      regMut = colSums(regMut)
      crtOut = myZtest(regMut, tExpD)    
      print(paste(crtOut,collapse=" "))
      outD = rbind(outD,crtOut)
    }
  }
}

outD = na.omit(outD)
outD[outD$V5 < 0.1,]
