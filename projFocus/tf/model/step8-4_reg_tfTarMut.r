rm(list=ls())

setRootd  = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysname']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/"
  }
  return(rootd)
}
rootd     = setRootd()

figd = paste(rootd, "/DATA/projFocus/report/Aug2014/fig",sep="")

source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

##func---

formatTestOut = function(test_aov){
  Res.DF_cnv_all = paste(test_aov$Res.Df, collapse="/")
  RSS_cnv_all = paste(round(test_aov$RSS,digits=2),collapse="/")  
  result = c(Res.DF_cnv_all= Res.DF_cnv_all, RSS_cnv_all= RSS_cnv_all, round(unlist(data.frame(test_aov)[2,3:6]),5))
  return(result)
}

anova_cnvVstfmut = function(inputdata){
  require(glmnet)
  tgene =   unlist(strsplit( tail(unlist(strsplit(x= inputdata, "/")),1), "\\."))[1]
  data = data.frame(t(read.table(inputdata, header=T, sep="\t",row.names=1,stringsAsFactors=F)))
  cnvall = data$cnv
  cnvSmp = which(cnvall !=0,arr.ind=T); cntCNVsmp = length(cnvSmp)
  cntTfmutSmp = length(which(rowSums(data[,-c(1,2)])!=0, arr.ind=T))
  
  if(cntTfmutSmp * 10 < cntCNVsmp){
    print(c(cntCNVsmp, cntTfmutSmp))
    p_prev = 1
    for(i in 1:100){
    data$cnv <- cnvall
    data$cnv[-sample(cnvSmp, cntTfmutSmp)] <- 0
  
#     data$cnv = sign(data$cnv)
    fitlm_selfcnv = lm(exp ~ cnv, data)
    fitlm_all = glmnet(exp ~ . ,alpha=1, data)
    tmp = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
#     print(res)
    if(as.numeric(tmp[6])<p_prev){
      res = tmp; p_rev = as.numeric(tmp[6])
    }
  }
  }else{
    
    fitlm_selfcnv = lm(exp ~ cnv, data)
    fitlm_all = glmnet(x= as.matrix(data[,-1]), y = data$exp ,alpha=1)
    
    res = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
  }
#   fitlm_muttf = lm(exp ~ ., data[,-2])
#   data$tfmut = ifelse(rowSums(data[,-c(1,2)]) > 0, 1, 0)
#   fitlm_allmerge = lm(exp ~ ., data)
#   data$exp = as.numeric(as.character(data[,1])) - as.numeric(as.character(fitlm_selfcnv$residuals))
#     
#   fitlm_all = cv.glmnet(y=data$exp, x = as.matrix(data[,-1]),nfolds=10)
#   plot(fitlm_all)

#     aovCNV = aov(exp ~ cnv, data)
#     aovAll = aov(exp ~., data)
#     print(aovCNV)
#     print(aovAll)
#     anova(aovCNV, aovAll)#   summary(aov(fitlm_tfmut))
#   summary(aov(fitlm_allmerge))
#   summary(aov(fitlm_all))

#   res = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
#   res = formatTestOut(anova(fitlm_selfcnv, fitlm_allmerge))
  return(c(tgene=tgene,cntCNVsmp, cntTfmutSmp,res))
}

anova_cnvVstfmut = function(inputdata){
  require(glmnet)
  tgene =   unlist(strsplit( tail(unlist(strsplit(x= inputdata, "/")),1), "\\."))[1]
  data = data.frame(t(read.table(inputdata, header=T, sep="\t",row.names=1,stringsAsFactors=F)))
#   data$cnv <- sign(data$cnv)
  fitlm_selfcnv = lm(exp ~ cnv, data)
  
  rss1 = sum(fitlm_selfcnv$residuals^2); 
  p1 = length(fitlm_selfcnv$coefficients ) - 1;
  
  fitlm_all = cv.glmnet(x= as.matrix(data[,-1]), y = data$exp ,alpha=1,nfolds=10, penalty.factor=c(0, rep(1,NCOL(data)-1)))
  # fitlm_all = cv.glmnet(x= as.matrix(data[,-1]), y = data$exp ,alpha=1,nfolds=10)
  
  modle2 = fitlm_all$glmnet.fit; mylambda = which(fitlm_all$lambda==fitlm_all$lambda.min)
  rss2 = sum( (data$exp- rowSums(modle2$beta[,mylambda] * X) - modle2$a0[mylambda] )^2)
  p2 = length(which(modle2$beta[,mylambda]!=0) )
  
  
  n = modle2$nobs
  F = ((rss1 - rss2)/(p2 -p1))/(rss2 /(n-p2))
  pval = ifelse(F>0,pf(F, p2-p1, n-p2, lower.tail=T), 1)

  model21 = formatTestOut(anova( lm(exp ~ ., data[,-2]), lm(exp~. , data)))
  print(c(tgene,F,pval, model21))
  return(c(tgene=tgene, F=F, pval=pval, model21 = model21 ) )
}
##func---end


###---main----
inputDir= "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/"
allInputFile =  vapply(list.files(path=inputDir,pattern="input"), FUN=function(x){
      paste(inputDir,"/",x,sep="")
},'aa')

result = sapply(allInputFile,anova_cnvVstfmut)

result = data.frame(t(result));rownames(result) = result[,1]; 
result = result[,-1]
length(which(as.numeric(as.character(result[,8]))<0.05))
result[(which(as.numeric(as.character(result[,8]))<0.05)),]
result['IFFO2',]
plot(as.numeric(as.character(result[,7])))

### test
inputdata="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/ERBB4.tgTFreg.input"
var.test(fitlm_selfcnv, fitlm_tfmut)
aov(fitlm_tfmut)
aov(fitlm_selfcnv)
result['HPS5',]


###====dignostic 
pdf(paste(figd,"/boxplot_cnvsmpVStfmutsmp_Aug09.pdf",sep=""))
boxplot(as.numeric(as.character(result$V2)), as.numeric(as.character(result$V3)), names=c("cnv_sample","tfmut_sample"), ylab = "event sample size", font=2,pch=19)
dev.off()

NROW(result[(which(as.numeric(as.character(result[,8]))<0.05)),])

kruskal.test(exp ~ cnv, data=data)
kruskal.test(exp~., data=data)
