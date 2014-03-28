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
usage = "Usage: Rscript grpLassoSNP.r --gene <gene name> --type <1/2/3>"
ERR = "ERROR here: "

#---funcs
trimSmpName = function (c) {
  return(vapply(c,FUN=function(x){substr(x,6,15)},'a'))
}

#---funce

##---init---
setRootd = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysnameß']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/projFocus/"
  }
  return(rootd)
}
rootd = setRootd()
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))


# #----getting command line parameters
# args = getArgs()
# plotflag = 0
# nperm    = 1000
# if(length(args) < 7 || is.null(args)){
#   print(paste(usage,example,sep="\n"))
#   print(args)
#   stop(paste(error,"wrong input parameter!"))
# }

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
wd          = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/CHEK1-temp/"
##---init
gene = tail(unlist(strsplit(wd,split="/|-",perl=T)),2)[1]
setwd(wd)
fexp  = "exp.mat"
fsmps = "samples.txt"
fsom  = "som.mat"
fcnv  = "cnv.mat"

##loading data
rawExp = read.delim2(fexp)
smpExp = colnames(rawExp)
colnames(rawExp) = smpExp
gExp = rownames(rawExp)

rawSmps = vapply(unlist(read.table("samples.txt")),as.character,'a')
smps = Reduce(intersect,list(rawSmps,smpExp))

dataExp = apply(subset(rawExp,select=smps),2,as.numeric  )
dataExp = t(dataExp[which(rowSums(dataExp)>0),])
dataExp = dataExp[apply(dataExp!=0,1,any),,drop=FALSE]
colnames(dataExp) = gExp


if (nrow(exp) <= 10) {
  require(car)
  library(MuMIn)
  require(MASS)
  require(glmnet)
  dataExp = as.data.frame(dataExp)
  #get candidate regulator for fix effect
  myFormula = as.formula(paste(colnames(dataExp)[1],"~", paste(colnames(dataExp)[-1],collapse="+"), sep=""))
  m = lm(myFormula, data =as.data.frame(dataExp))
  fm  = stepAIC(m,trace=T)
  sfm = summary(fm)
  
  ## select candidate driver
  require(nlme)
  dataMat = dataExp[1:5,1:5]
  mat2nlmeData = function(dataMat){
  require(reshape)
  dataMat$sbj = rownames(dataMat)
    dout <- reshape(dataMat, 
                 varying = colnames(dataMat[-c(1,ncol(dataMat))]), 
                 v.names = "regulator.exp",
                 new.row.names = paste(rownames(dataMat),1:25,sep="."),
                 direction = "long")
    return(dout)    
  }
  
  fitCoeff = as.data.frame(sfm$coefficients[-grep("Intercept",rownames(sfm$coefficients)),])
  colnames(fitCoeff) = c(colnames(fitCoeff)[1:3],"p.value")
  condsCoeff = fitCoeff[which(fitCoeff$p.value < 0.05),]
  cands = rownames(condsCoeff[order(condsCoeff[,4]),])
  
  myFix = as.formula(paste(colnames(dataExp)[1], "~",
                               paste(cands,collapse="*"),sep=""))
  myRandom = as.formula(paste("~ ",paste(colnames(dataExp)[-1],collapse="+"),sep=""))
  myCluster= rownames(dataExp)
  exp.lme = lme(fixed=myFix,random=myRandom,data=dataExp)
  for (cand in seq_along(cands)){
    cand = cands[4]
  }
  
  
  
}
else{
  #feature selection
  exp = t(dataExp[,-1])
  doGrouping(){
    ##svd transform, take 80% variance, do kmeans, and then output grouping for each variable
        vars = rownames(exp)
        dexp = dist(exp)
        group         = getGroup(exp,dexp)
        names(group)  = vars
    }
  
  cntsample = nrow(dataExp) ; cntReg = ncol(dataExp) -1;
  require(grpreg)
  myexp.gl = regfit(as.matrix(dataExp[,-1]),as.matrix(dataExp[,1]),group)
  
  print("Doing regression...")
  nperm     = 1000
  fitpermut = regfitPermu(dataExp, group, 1000, plotflag=0) 
  coeff = myexp.gl$beta[order(abs(myexp.gl$beta),decreasing=T)]
  coeff = coeff[-grep("Intercept",names(coeff))]
  coeff = sort(coeff[coeff!=0],decreasing=T)
  cddts = names(coeff)
  
}

###-------


fixedFactor = names(sFit$coefficients[-1])
library(nlme)


require("ICC")
require(fpc)
for (i in seq(min(dexp) * 1.5, max(dexp) * 0.8,by=min(dexp) * 1.5)){
  fitDB = dbscan(data=dexp,eps=15,showplot=2,method="raw")
  plot(fitDB,exp)
  fitDB$cluster
  cor.test(exp[rownames(exp)[fitDB$cluster==0],] (dataExp[,1]))
  t.test(exp[rownames(exp)[fitDB$cluster==1],],(dataExp[,1]))
  icc = ICCbare(x=fitDB$cluster, y = rownames(exp),data = t(exp))  
  table(fitDB$cluster)
  cor()
  print()
}

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
     par(mar=c(0,0,0,0))
     heatmap.2(dataExp,col=mycol,trace='none')
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
