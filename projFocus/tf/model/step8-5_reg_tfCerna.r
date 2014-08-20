##status: under work

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

###----func-------
formatTestOut = function(test_aov){
  Res.DF_cnv_all = paste(test_aov$Res.Df, collapse="/")
  RSS_cnv_all = paste(round(test_aov$RSS,digits=2),collapse="/")  
  result = c(Res.DF_cnv_all= Res.DF_cnv_all, RSS_cnv_all= RSS_cnv_all, round(unlist(data.frame(test_aov)[2,3:6]),5))
  return(result)
}

orderSampleLable = function(dataExp){
  dataExp$sample <- as.character(dataExp$sample)
  #Then turn it back into an ordered factor
  dataExp$sample <- factor(dataExp$sample, levels=unique(dataExp$sample))
  return(dataExp)
}

noseHeatmap = function(inputfile){
#   inputfile = paste(inputDir, "/input_", tgene, sep="")
  tgene = tail(unlist(strsplit(inputfile,"_")),1)
  data = read.csv2(inputfile, sep="\t", header=T,stringsAsFactors=F)
  data = data[order(as.numeric(as.character(data$cTarExp))),]
  allSmp = data$X; rownames(data) <- allSmp
  data = data[,colSums(data!=0)!=0]
  
  allFeatures = colnames(data)
  data$sample = rownames(data)
  ftr.cRegExp = allFeatures[grep("cRegExp", allFeatures)]
  ftr.cRegCNV = allFeatures[grep("cRegCNV", allFeatures)]
  ftr.cRegTFmut = allFeatures[grep("cRegTFmut", allFeatures)]
  data[,c('cTarExp', ftr.cRegExp, ftr.cRegCNV, ftr.cRegTFmut)] = apply(data[,c('cTarExp', ftr.cRegExp, ftr.cRegCNV, ftr.cRegTFmut)],2,as.numeric)
  
  library(reshape); library(ggplot2) ;library(scales);library(RCurl);  library (grid)
  
  ####---underDevelopment
  
  ### ceRNA target expression
  data.cTarExp = data[,c('sample', 'cTarExp')]
  data.cTarExp$variable = rep("   cTarExp ",NROW(data)); colnames(data.cTarExp) = c("sample", 'value', 'variable')
  data.cTarExp$value <- as.numeric(as.character(data.cTarExp$value))
  data.cTarExp <- orderSampleLable(data.cTarExp)
  data.cTarExp$rescale <- rescale(data.cTarExp$value,to=c(-1,1))
  head(data.cTarExp)
  
  ### ceRNA regulator expression
  data.cRegExp = data[,c('sample', ftr.cRegExp)]
  if (length(ftr.cRegExp) > 1){
    data.cRegExp[, ftr.cRegExp] = apply(data.cRegExp[,ftr.cRegExp],2, function(x){rescale(x,to=c(-1,1))})
  }else if (length(ftr.cRegExp) == 1){
    data.cRegExp[, ftr.cRegExp] = rescale( data.cRegExp[, ftr.cRegExp], to = c(-1,1))
  }else {
    print("no cRegExpression information!")
    return(0)
  }
  data.cRegExp.m <- melt(data.cRegExp, id.vars='sample')
  head(data.cRegExp.m) 
  data.cRegExp.m <- orderSampleLable(data.cRegExp.m)
  head(data.cRegExp.m)

  
  ## ceRNA regulator's CNV 
  data.cRegCNV = data[,c('sample', ftr.cRegCNV)]
  data.cRegCNV = melt(data.cRegCNV, id.vars='sample')
  data.cRegCNV$value <- factor(sign(as.numeric(as.character(data.cRegCNV$value))), levels = c(-1,0,1))
  data.cRegCNV <- orderSampleLable(data.cRegCNV)
  head(data.cRegCNV)
  
  ## ceRNA regulator's TF mutation profile (collapsed)
  data.cRegTFmut  = data[,c('sample',ftr.cRegTFmut)]; data.cRegTFmut[,ftr.cRegTFmut] <- sign(data.cRegTFmut[,ftr.cRegTFmut])
  data.cRegTFmut.m = melt(data.cRegTFmut,id.vars='sample')
  data.cRegTFmut.m$value <- factor(data.cRegTFmut.m$value, levels=c(0,1))
  data.cRegTFmut.m <- orderSampleLable(data.cRegTFmut.m)
  
  head(data.cRegTFmut.m)
  
  ### ploting
    base_size <- 6
    p <- ggplot(data.cTarExp, aes(variable, sample) )+ geom_tile(aes(fill = rescale), color = 'white') 
    p1 <- p + scale_fill_gradient2(low='green',mid='white',high='orange') + 
          theme_grey(base_size = base_size) + theme(axis.title=element_blank(), 
              axis.text.y = element_blank(),
#               axis.text.x = element_blank(),            
              legend.position="none",
#               axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
              axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=16))


    p <- ggplot(data.cRegExp.m, aes(variable, sample)) +   geom_tile(aes(fill = value), color = 'white')
    p2 <- p + scale_fill_gradient2(low='green',mid='white',high='orange') + 
            theme_grey(base_size = base_size) +  theme(axis.title=element_blank(), 
            axis.text.y = element_blank(),
#             axis.text.x = element_blank(), 
            legend.position="none" ,
#             axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=6), axis.text.x  = element_text(angle=90, vjust=0.5, size=6))
        
    p <- ggplot(data.cRegCNV, aes(variable, sample)) +   geom_tile(aes(fill = value), color = 'white')
    p3 <- p + scale_fill_manual(values=c("lightblue","white",colors()[119])) +  
#           theme_grey(base_size = base_size) +   
            theme(axis.title=element_blank(), 
                  axis.text.y = element_blank(),
                    legend.position="none", 
#                     axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=6), axis.text.x  = element_text(angle=90, vjust=0.5, size=6))
      
    p <- ggplot(data.cRegTFmut.m, aes(variable, sample)) +  geom_tile(aes(fill = value), colour =   "white") 
    p4 <- p + scale_fill_manual(values=c('white',"red")) + 
            theme(axis.title=element_blank(),  
                  axis.text.y = element_blank(), 
                    legend.position="none",
#             axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=6), axis.text.x  = element_text(angle=90, vjust=0.5, size=6))
    
     
  len1 = length(unique(data.cTarExp$variable)); len2 = length(unique(data.cRegExp.m$variable)); len3= length(unique(data.cRegCNV$variable)); len4=length(unique(data.cRegTFmut.m$variable)); 
  lensum = len2 + len3 + len4
  pdf(paste(figd,"/heatmap_diagnois_regCerna_Aug17_",tgene, ".pdf",sep=""),width=15,height=18)
  layt<-grid.layout(nrow=1,ncol=4,widths=c(1/20,len2/lensum * 19/20,len3/lensum * 19/20,len4/lensum * 19/20))
  
  layt<-grid.layout(nrow=1,ncol=4,widths=c(1/15,3/10,3/10,3/10))
  
  #   grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p1,  vp = viewport(layout.pos.row=1,layout.pos.col=1))
  print(p2, vp = viewport(layout.pos.row=1,layout.pos.col=2))
  print(p3, vp = viewport(layout.pos.row=1,layout.pos.col=3))
  print(p4, vp = viewport(layout.pos.row=1,layout.pos.col=4))
  dev.off()

}

svReg = function(inputfile) {
  tgene = tail(unlist(strsplit(inputfile,"_")),1)
  
  data = read.csv2(inputfile, sep="\t", header=T,stringsAsFactors=F)
  data = data[order(as.numeric(as.character(data$cTarExp))),]
  allSmp = data$X; rownames(data) <- allSmp
  data = data[,colSums(data!=0)!=0]
  
  
  allFeatures = colnames(data)
  ftr.cTarExp = allFeatures[grep("cTarExp", allFeatures)]
  ftr.cRegExp = allFeatures[grep("cRegExp", allFeatures)]
  ftr.cRegCNV = allFeatures[grep("cRegCNV", allFeatures)]
  ftr.cRegTFmut = allFeatures[grep("cRegTFmut", allFeatures)]
  
  if (length(ftr.cRegTFmut) == 0 ){
    return(c(tgene, length(ftr.cRegCNV), 0, NROW(data), NA, NA))
  }
  
  data[,c(ftr.cRegTFmut, ftr.cRegCNV)] <- sign(apply(data[,c(ftr.cRegTFmut, ftr.cRegCNV)],2, as.numeric))
  data[,c(ftr.cTarExp, ftr.cRegExp)] <- apply(data[,c(ftr.cTarExp, ftr.cRegExp)], 2, as.numeric)
  
  
  
  ### plot
  # noseHeatmap(data)
  
  data.mod1 = as.matrix(data[,c(ftr.cTarExp, ftr.cRegCNV)])
  data.mod2 = as.matrix(data[,c(ftr.cTarExp, ftr.cRegCNV, ftr.cRegTFmut)])
  
  require(e1071)
  cost1 = 100; cost2 = 1000; gamma = 0.001
  mod1 = svm(cTarExp ~ ., data = data.frame(data.mod1), type="eps-regression",  kernel = "polynomial", cost = cost1, gamma = gamma, cross=10)
  mod2 = svm(cTarExp ~ ., data = data.frame(data.mod2), type="eps-regression",  kernel = "polynomial", cost = cost2, gamma = gamma, cross = 10)
#   testReg1 = tune(svm, cTarExp ~ ., data = data.frame(data.mod1), ranges=list(type="eps-regression",  kernel = "polynomial", cost = 2^(5:10), gamma = 10^(-5:1)),
#                       tunecontrol = tune.control(sampling = "fix")  )
#   
#   
#   testReg2 = tune(svm, cTarExp ~ ., data = data.frame(data.mod2), ranges=list(type="eps-regression",  kernel = "polynomial", cost = 2^(5:10), gamma = 10^(-5:1)),
#                       tunecontrol = tune.control(sampling = "fix"))
#   print(c(testReg1$best.parameter, testReg2$best.parameter))
#   mod1 = testReg1$best.model
#   mod2 = testReg2$best.model
  
  rss1 = sum(mod1$residuals)^2 ; rss2 = sum(mod2$residuals)^2
  p1 = NCOL(data.mod1) ; p2 = NCOL(data.mod2)
  n = NROW(data)
  
  ##alternative hypo: mod2 is bettern than mod1, rss1 is smaller than rss2
  f = ((rss1-rss2)/(p2-p1))/(rss2/(n-p2))
  f.pval = ifelse(f>0,pf(f, p2-p1, n-p2, lower.tail=T), 1)
  out = c(tgene, p1, p2-p1, n, f, f.pval)
  print(out)
  return(out)
}

###----main-------
require(scales)
inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/data" , sep="")
inputfileLst = paste(inputDir, "/", list.files(inputDir), sep="")
output  = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/result_", sep="")

tgene = 'EN2'
inputfile = paste(inputDir, '/', 'input_', tgene, sep="")
# inputfile = sample(inputfileLst, 1)
# = function(inputfile) {
#   tgene = tail(unlist(strsplit(inputfile,"_")),1)
  
  data = read.csv2(inputfile, sep="\t", header=T,stringsAsFactors=F)
  data = data[order(as.numeric(as.character(data$cTarExp))),]
  allSmp = data$X; rownames(data) <- allSmp
  data = data[,colSums(data!=0)!=0]
  
  
  allFeatures = colnames(data)
  ftr.cTarExp = allFeatures[grep("cTarExp", allFeatures)]
  ftr.cRegExp = allFeatures[grep("cRegExp", allFeatures)]
  ftr.cRegCNV = allFeatures[grep("cRegCNV", allFeatures)]
  ftr.cRegTFmut = allFeatures[grep("cRegTFmut", allFeatures)]
  
  if (length(ftr.cRegTFmut) == 0 ){
    return(c(tgene, length(ftr.cRegCNV), 0, NROW(data), NA, NA))
  }
  
  data[,c(ftr.cRegTFmut, ftr.cRegCNV)] <- sign(apply(data[,c(ftr.cRegTFmut, ftr.cRegCNV)],2, as.numeric))
  data[,c(ftr.cTarExp, ftr.cRegExp)] <- apply(data[,c(ftr.cTarExp, ftr.cRegExp)], 2, as.numeric)
  
  
  
  ### plot
  # noseHeatmap(data)
  
  data.mod1 = data[,c(ftr.cTarExp, ftr.cRegCNV)]
  data.mod2 = data[,c(ftr.cTarExp, ftr.cRegCNV, ftr.cRegTFmut)]
  kruskal.test(cTarExp ~ ., data = data.mod1)
  kruskal.test(cTarExp ~ ., data = data.mod2)

  require(e1071)
  cost1 = 100; cost2 = 1000; gamma = 0.001
  mod1 = svm(cTarExp ~ ., data = data.frame(data.mod1), type="eps-regression",  kernel = "polynomial", cost = cost1, gamma = gamma, cross=10)
  mod2 = svm(cTarExp ~ ., data = data.frame(data.mod2), type="eps-regression",  kernel = "polynomial", cost = cost2, gamma = gamma, cross = 10)
  #   testReg1 = tune(svm, cTarExp ~ ., data = data.frame(data.mod1), ranges=list(type="eps-regression",  kernel = "polynomial", cost = 2^(5:10), gamma = 10^(-5:1)),
  #                       tunecontrol = tune.control(sampling = "fix")  )
  #   
  #   
  #   testReg2 = tune(svm, cTarExp ~ ., data = data.frame(data.mod2), ranges=list(type="eps-regression",  kernel = "polynomial", cost = 2^(5:10), gamma = 10^(-5:1)),
  #                       tunecontrol = tune.control(sampling = "fix"))
  #   print(c(testReg1$best.parameter, testReg2$best.parameter))
  #   mod1 = testReg1$best.model
  #   mod2 = testReg2$best.model
  
  rss1 = sum(mod1$residuals)^2 ; rss2 = sum(mod2$residuals)^2
  p1 = NCOL(data.mod1) ; p2 = NCOL(data.mod2)
  n = NROW(data)
  
  ##alternative hypo: mod2 is bettern than mod1, rss1 is smaller than rss2
  f = ((rss1-rss2)/(p2-p1))/(rss2/(n-p2))
  f.pval = ifelse(f>0,pf(f, p2-p1, n-p2, lower.tail=T), 1)
  out = c(tgene, p1, p2-p1, n, f, f.pval)
  print(out)
  return(out)
# }



# ---------------



resDF = sapply(inputfileLst,svReg)

resDF = t(resDF)
write.table(resDF, file= paste(output, "_Aug182014",sep=""),col.names=F, quote=F, sep="\t" )
save.image(paste(output, ".Rda", sep=""))


pcut = 0.01
plot(which(as.numeric(resDF[,6]) <= pcut))

colnames(resDF) <- c('cTar', 'cnt_cRegCNV', 'cnt_cRegTFmut', 'cnt_smp','F', 'F.pvalue')
resDF = data.frame(resDF)
row.names(resDF) <- resDF$cTar
resDF[,c('cnt_cRegCNV', 'cnt_cRegTFmut', 'cnt_smp','F', 'F.pvalue')] <- apply(resDF[,c('cnt_cRegCNV', 'cnt_cRegTFmut', 'cnt_smp','F', 'F.pvalue')], 2, function(x){as.numeric(as.character(x))})


pdf(paste(figd,"/scatterplot_regCernaTF_SVR_Aug19.pdf", sep=""))
plot.y = resDF$F.pvalue
plot.y[which(plot.y == 0)] <- 1e-300
plot(x = resDF$cnt_cRegTFmut/resDF$cnt_cRegCNV, y = -log10(plot.y), col = 'blue', 
     xlab = "cnt_mut/cnt_cnv", ylab = "-log10(p-value)", main ="P value of F test\n(SVR)" , 
     font = 2)

# plot(x = resDF$cnt_cRegTFmut, y = -log10(resDF$F.pvalue), col = 'blue', 
#      xlab = "cnt_mut", ylab = "-log10(p-value)", main ="P value of F test\n(SVR)" , 
#      font = 2)

abline(h=2,col="red", lwd = 4);
text(x=40, y=300, labels=paste("Significant:  ", round(length(which(resDF$F.pvalue<=pcut))/NROW(resDF) * 100,digits=3), "%\n       ", 
                              length(which(resDF$F.pvalue<=0.01)), " / ", NROW(resDF), sep=""), font = 2)
dev.off()
 
resDF.sig = resDF[which(resDF[,6]<=pcut),]
resDF.nosig = resDF[which(resDF[,6]>pcut),]


### analysis not significant ones
# targs = inputfileLst[3:20]
# for (inputfile in targs ) {
#     noseHeatmap(inputfile)
# }

cTars.fail = as.character(unlist(resDF.nosig$cTar))

figd = paste(rootd, "/DATA/projFocus/report/Aug2014/fig/regCernaFail/",sep="")

for (cTar in sample(cTars.fail,size=max(length(cTars.fail),20))){
  inputfile= paste(inputDir, "/input_", cTar, sep="")
  print(c(cTar, inputfile))
  noseHeatmap(inputfile)
}

cTars.succ = as.character(unlist(resDF.sig$cTar))
figd = paste(rootd, "/DATA/projFocus/report/Aug2014/fig/regCernaSucc/",sep="")
plot.cTars = sample(cTars.succ,size=10,replace=F)

for (cTars.crt in plot.cTars){
  inputfile= paste(inputDir, "/input_", cTars.crt, sep="")
  print(c(cTars.crt, inputfile))
   noseHeatmap(inputfile)
}

resDF.sig[order(resDF.sig$F.pvalue),][1000:1030,]
top.cTars = row.names(resDF.sig[order(resDF.sig$F.pvalue),][1000:1010,])
top.cTars =c("STAT3", "TBX3")
for (cTars.crt in top.cTars){
  inputfile= paste(inputDir, "/input_", cTars.crt, sep="")
  print(c(cTars.crt, inputfile))
  noseHeatmap(inputfile)
}
