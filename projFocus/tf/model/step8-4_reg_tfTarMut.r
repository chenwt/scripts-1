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
#   data$cnv = sign(data$cnv)
  fitlm_selfcnv = lm(exp ~ cnv, data)
  fitlm_all = lm( exp ~ . , data)
#   fitlm_muttf = lm(exp ~ ., data[,-2])
#   data$tfmut = ifelse(rowSums(data[,-c(1,2)]) > 0, 1, 0)
#   fitlm_allmerge = lm(exp ~ ., data)
#   data$exp = as.numeric(as.character(data[,1])) - as.numeric(as.character(fitlm_selfcnv$residuals))
#     
#   fitlm_all = cv.glmnet(y=data$exp, x = as.matrix(data[,-1]),nfolds=10)
#   plot(fitlm_all)

    aovCNV = aov(exp ~ cnv, data)
    aovAll = aov(exp ~., data)
    print(aovCNV)
    print(aovAll)
    anova(aovCNV, aovAll)#   summary(aov(fitlm_tfmut))
#   summary(aov(fitlm_allmerge))
#   summary(aov(fitlm_all))

  res = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
#   res = formatTestOut(anova(fitlm_selfcnv, fitlm_allmerge))

  return(c(tgene=tgene,res))
}

rss1 = sum(fitlm_selfcnv$residuals^2); 
p1 = length(fitlm_selfcnv$coefficients ) - 1;

rss2 = sum(fitlm_all$residuals^2)
p2 = length(fitlm_all$coefficients) - 1

sum(fitlm_tfmut$residuals^2)
length(fitlm_muttf$coefficients) - 1

n = fitlm_all$df.residual
((rss1 - rss2)/(p2 -p1))/(rss2 /(n-p2))
##func---end

inputDir= "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/"
allInputFile =  vapply(list.files(path=inputDir,pattern="input"), FUN=function(x){
      paste(inputDir,"/",x,sep="")
},'aa')

result = sapply(allInputFile,anova_cnvVstfmut)
result = data.frame(t(result));rownames(result) = result[,1]; 
result = result[,-1]
length(which(as.numeric(as.character(result[,6]))<0.05))
result['IFFO2',]
plot(as.numeric(as.character(result[,7])))
### test
inputdata="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/IFFO2.tgTFreg.input"
var.test(fitlm_selfcnv, fitlm_tfmut)
aov(fitlm_tfmut)
aov(fitlm_selfcnv)
result['HPS5',]
