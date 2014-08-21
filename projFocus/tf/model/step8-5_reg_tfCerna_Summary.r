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
figd = paste(rootd, "/DATA/projFocus/report/Aug2014/fig/regBoot",sep="")
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))


#----getting command line parameters
#--getting arguments
# args = getArgs()
# 
# if(length(args) < 1 || is.null(args)){
#     print(paste(usage,example,sep="\n"))
#   print(args)
#     stop(paste(error,"wrong input parameter!"))
# }
# 
# setwd(system("pwd",intern=T))
# tgene    = args[['ctar']]
# # imputfile = 
# nboot	 = as.integer(args[['nboot']])
# # output   = args[['out']]

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

noseHeatmap = function(tgene){
  inputfile = paste(inputDir, "/input_", tgene, sep="")
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
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
        
    p <- ggplot(data.cRegCNV, aes(variable, sample)) +   geom_tile(aes(fill = value), color = 'white')
    p3 <- p + scale_fill_manual(values=c("lightblue","white",colors()[119])) +  
#           theme_grey(base_size = base_size) +   
            theme(axis.title=element_blank(), 
                  axis.text.y = element_blank(),
                    legend.position="none", 
#                     axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
      
    p <- ggplot(data.cRegTFmut.m, aes(variable, sample)) +  geom_tile(aes(fill = value), colour =   "white") 
    p4 <- p + scale_fill_manual(values=c('white',"red")) + 
            theme(axis.title=element_blank(),  
                  axis.text.y = element_blank(), 
                    legend.position="none",
#             axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
    
     
  len1 = length(unique(data.cTarExp$variable)); len2 = length(unique(data.cRegExp.m$variable)); len3= length(unique(data.cRegCNV$variable)); len4=length(unique(data.cRegTFmut.m$variable)); 
  lensum = len2 + len3 + len4
  pdf(paste(figd,"/heatmap_nose_regCerna_Aug21_",tgene, ".pdf",sep=""),width=15,height=18)
  layt<-grid.layout(nrow=1,ncol=4,widths=c(1/20,len2/lensum * 19/20,len3/lensum * 19/20,len4/lensum * 19/20))
  
  layt<-grid.layout(nrow=1,ncol=3,widths=c(1/21,len3/lensum * 10/21,len4/lensum * 10/21))

  #   grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p1,  vp = viewport(layout.pos.row=1,layout.pos.col=1))
#   print(p2, vp = viewport(layout.pos.row=1,layout.pos.col=2))
  print(p3, vp = viewport(layout.pos.row=1,layout.pos.col=2))
  print(p4, vp = viewport(layout.pos.row=1,layout.pos.col=3))
  dev.off()

}

calKS = function(inputfile){
  tgene = tail(unlist(strsplit(inputfile,"_")),1)
  if( ! file.exists(inputfile) | file.info(inputfile)$size == 0 ){
    return(rep(NA,5))
  }
  print(inputfile)
  data = read.csv(inputfile, stringsAsFactors= F, header=F,sep="\t")
  colnames(data ) = c("btname","tgene", "mod1.mse", "mod1.r2",
                      "mod2.mse", "mod2.r2",
                      "rss1", "rss2",
                      "cntCNV", "cntMut", "cntSmp", "F", "FPval")
  
  ks.mse = ks.test(as.numeric(data$mod1.mse), as.numeric(data$mod2.mse),alternative='less')
  ks.r2 = ks.test(as.numeric(data$mod1.r2), as.numeric(data$mod2.r2),alternative='greater')
  out= sprintf("%s\t%5.3f\t%5.3e\t%5.3f\t%5.3e" ,tgene,ks.mse$statistic,ks.mse$p.value, ks.r2$statistic, ks.r2$p.value)
  return(out)
}

#### main------
require(scales)
inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/runAug21/" , sep="")

output  = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/summary_runAug21", sep="")


inputfileLst = paste(inputDir, "/", list.files(inputDir,pattern="result"), sep="")
# inputfile = paste(inputDir, '/', 'result_', tgene, sep="")
resDF = sapply(inputfileLst, calKS)

resDF = t(resDF)
names(resDF) <- "cTarKS.mse\tKS.mse.pval\tKS.r2\tKS.r2.pval"
write(resDF, file=output )

resDF = read.table(output, sep="\t",stringsAsFactors=F)
pcut = 0.001

NROW(resDF[which(as.numeric(resDF[,3]) <= pcut),])/ NROW(resDF)
resDF.sig = resDF[which(as.numeric(resDF[,3]) <= pcut),]
resDF.sig = resDF.sig[order(resDF.sig$V3),]
head(resDF.sig)
###----plot density
inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/runAug21/" , sep="")
plotMSEDensity = function(tgene){
#   tgene = tail(unlist(strsplit(inputfile,"_")),1)
  inputfile = paste(inputDir, '/', 'result_', tgene, sep="")
  
  if( ! file.exists(inputfile) | file.info(inputfile)$size == 0 ){
    return(rep(NA,5))
  }
  print(inputfile)
  data = read.csv(inputfile, stringsAsFactors= F, header=F,sep="\t")
  colnames(data ) = c("btname","tgene", "mod1.mse", "mod1.r2",
                      "mod2.mse", "mod2.r2",
                      "rss1", "rss2",
                      "cntCNV", "cntMut", "cntSmp", "F", "FPval")
  data[,c(3,4,5,6)] = apply(data[,c(3,4,5,6)], 2, as.numeric)
  
  pdf(paste(figd,"/density_",tgene, ".pdf", sep=""))
  
    plot(0, type="n", ylab = "density", xlab="",xlim = c(min(data$mod1.mse, data$mod2.mse),max(data$mod1.mse, data$mod2.mse)),
         ylim=c(0,max(density(as.numeric(data$mod1.mse))$y, density(as.numeric(data$mod2.mse))$y)),
         main = paste(tgene, "MSE"))
    lines(density(as.numeric(data$mod1.mse)),col='blue', lwd = 2,)
    lines(density(as.numeric(data$mod2.mse)), col = 'red', lwd=2)
    legend('topright',legend=c("model1 MSE", "model2 MSE"),lty=1, col=c("blue","red"), bty ='n', cex = .75, lwd=2)
  dev.off()
}
plotMSEDensity("AFF3")
plotMSEDensity("ABCA1")
plotMSEDensity("CSF1")


inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/data/" , sep="")
noseHeatmap('CSF1')
noseHeatmap('AFF3')

resDF.nosig = resDF[which(as.numeric(resDF[,3])  > pcut),]
resDF.nosig = resDF.nosig[order(resDF.nosig$V3,decreasing=T),]
head(resDF.nosig)
