##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <string: regulator gene name > 
#output: <file: > 
#Usage: Rscript $scriptName 
#Description: require package irr, grpreg,gplots
#             Major change: analysis for one regulator 
#                           formalized model                           
#TODO: develop plot option!
example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript 
    /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/grpLassoSNP.r
    --gene SMAD4 
    --type 1 
    --wd /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/temp-ESR1 "
usage = "Usage: Rscript step2-3_regulatorGrpLasso.R --gene <gene name> --type <1/2/3>"
ERR = "ERROR here: "

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  rootd= "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
 }else if(sysInfo['sysname']=="Linux" ){
  rootd= "/ifs/home/c2b2/ac_lab/jh3283/"
}
source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/myR/general.r")
source(jxy(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R"))
source(jxy(rootd, "/scripts/projFocus/ceRNA/model/grpLassoFunctions.R"))
setwd(jxy(rootd,"/DATA/projFocus/result/02022014/model/"))

#----getting command line parameters
# args = getArgs()
# plotflag = 0
# nperm    = 1000
# if(length(args) < 7 || is.null(args)){
#   print(paste(usage,example,sep="\n"))
#   print(args)
#   stop(paste(error,"wrong input parameter!"))
# }
# 
# setwd(system("pwd",intern=T))
# cwd         = getwd()
# gene        = args['gene']
# type        = args['type']
# wd          = args['wd']
# 
# plotflag    = as.integer(args['plot']) #optional 
# nperm       = as.integer(args['nperm'])
# output      = paste(cwd,"/",args['out'], sep="")
# print(paste("current working directory:",cwd))

#--test
cwd         = getwd()
wd          = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/temp-ESR1"
gene        = "TCF4"
out         = "grplasso"
output      = paste(wd,"/",out, sep="")

##-------------prepareData
setwd(wd)
inputsnp    = paste(wd,"/",gene, ".snp",sep="")
inputexp    = paste(wd,"/",gene, ".exp",sep="")
inputcnv    = paste(wd,"/",gene, ".cnv",sep="")
inputsom    = paste(wd,"/", ".som",sep="")
inputSample = paste(wd,"/", gene, ".sample",sep="") ##<<====HERE
inputRegulator = paste(wd,"/regulator",sep="")

##---testleve---
smps = unlist(strsplit(unlist(read.table(inputSample,header=T,stringsAsFactors=F)[1,2]),";"))
dataSom = t(as.data.frame(runif(length(smps),0,1)))
colnames(dataSom) = smps
rownames(dataSom) = "ESR1:1:111:112"
#---

##----------------------------##----------------------------
##load exp data
dataExp = formatData(inputFile=inputexp,t='exp')
dataExp = normalize(dataExp)
outtxt = paste(output,gene,".txt",sep="")
outplot = paste(output,gene,".pdf",sep="")

if (plotflag == 1){
     print("plotting...")
     pdf(outplot)
     require(gplots)
     mycol = bluered(256)
     plot(density(dataExp),main = "Gene Expression Density")
     blankPlot(100,1)
     image(dataExp,col=mycol)
     axis(2, at=unlist(sapply(seq(0.5,92.5),FUN=function(x){c(x,0))})), labels=rownames(dataExp), lty=, col=, las=2) 
     
    }
dataSnp = formatData(inputFile=inputsnp,t='snp')

if (dataCnt$snp > 1){
  print("grouping snps...")
  kcDist          = kappaDist(dataSnp) 
  kcDist_svd      = getSVD(kcDist)
  group           = getGroup(kcDist_svd,kcDist)
  if(plotflag == 1){
    image(z=kcDist,col=mycol,main="image: KC similarity SNP")
    heatmap.2(kcDist,trace="none",col=mycol,dendrogram="none", main="heatmap kc similarity")
    image(z=kcDist_svd,col=mycol,main="image: KC similarity SNP after svd")
    row.names(kcDist_svd) = rownames(kcDist)
    heatmap.2(kcDist_svd,col=mycol,trace='none',dendrogram="none", main="heatmap: KC similarity SNP after svd",labCol=" ")
  }
}else{
  print("only one snp...")
  group   = as.integer(1)
  names(group) = colnames(dataSnp)
}  
switch(type,'1'={
    ##model 1 exp = snp + som + cnv
    print("model exp~snp+som+cnv")
    dataCnv = formatData(inputcnv,t='cnv')
    dataSom = formatData(inputsom,t='som')
    cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntcnv = ncol(dataCnv); cntsom = ncol(dataSom)
    group     = c(som = rep(1,cntsom), cnv = rep(2,cntcnv), group + 2)
    
    data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntcnv+cntsom))
    colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataSom),colnames(dataCnv))
    row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames
    
    data_merge[,1] = dataExp
    data_merge[,2:(1+cntsnp)] = dataSnp
    data_merge[,(1+cntsnp+1):(1+cntsnp+cntsom)] = dataSom
    data_merge[,(1+cntsnp+cntsom+1):(1+cntsnp+cntsom+cntcnv)] = dataCnv

    },'2'={    
      ##model 2 exp = snp + soｍ 
      print("model exp~snp + som")
      dataSom = formatData(inputsom,t='som')
      cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntsom = ncol(dataSom)
      group     = c(som = rep(1,cntsom), group + 1)
      
      data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntsom))
      colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataSom))
      row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames
      
      data_merge[,1] = dataExp
      data_merge[,2:(1+cntsnp)] = dataSnp
      data_merge[,(1+cntsnp+1):(1+cntsnp+cntsom)] = dataSom
        
    },'3'={
     ## model 3 exp = snp + cnv 
      print("model exp~snp + cnv")
      dataCnv = formatData(inputcnv,t='cnv')
      cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp); cntcnv = ncol(dataCnv)
      group     = c( cnv = rep(1,cntcnv), group + 1)
      
      data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp+cntcnv))
      colnames(data_merge) = c('exp',colnames(dataSnp),colnames(dataCnv))
      row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames

      data_merge[,1] = dataExp
      data_merge[,2:(1+cntsnp)] = dataSnp
      data_merge[,(1+cntsnp+1):(1+cntsnp+cntcnv)] = dataCnv      
    }, '4'={
      print("model exp ~ snp")
      cntsample = nrow(dataExp) ; cntsnp = ncol(dataSnp);
      data_merge = as.data.frame(matrix(NA,nrow=cntsample,ncol=1+cntsnp))
      colnames(data_merge) = c('exp',colnames(dataSnp))
      row.names(data_merge) = row.names(dataExp) ## assume all data matrix with same rownames    
      data_merge[,1] = dataExp
      data_merge[,2:(1+cntsnp)] = dataSnp
    })
  
require(grpreg)
print("Doing regression...")
nperm     = nperm
fitpermut = regfitPermu(data_merge, group, nperm, plotflag=plotflag) 

##output
print("writing output...")
write.table(t(as.matrix(c(paste("gene","RSS","npermu","pvalue",sep=":"),paste(genename,fitpermut$RSS, nperm, fitpermut$pvalue,sep=":")))),
            outtxt,
            quote=F,col.names=F,sep="\t",row.names = F)
write.table(as.matrix(sort(fitpermut$beta[-1])),
            outtxt,
            append=T,
            col.names=F, quote=F,sep="\t")
print(paste("#-----Done",genename))
