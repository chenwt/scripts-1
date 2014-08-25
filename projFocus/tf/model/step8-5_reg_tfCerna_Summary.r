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
  ftr.cTarExp = allFeatures[grep("cTarExp", allFeatures)]
  ftr.cTarCNV = allFeatures[grep("cTarCNV", allFeatures)]
  ftr.cRegExp = allFeatures[grep("cRegExp", allFeatures)]
  ftr.cRegCNV = allFeatures[grep("cRegCNV", allFeatures)]
  ftr.cRegTFmut = allFeatures[grep("cRegTFmut", allFeatures)]
  data[,c(ftr.cTarExp, ftr.cRegExp, ftr.cTarCNV, ftr.cRegCNV, ftr.cRegTFmut)] = 
        apply(data[,c(ftr.cTarExp, ftr.cRegExp, ftr.cTarCNV,ftr.cRegCNV, ftr.cRegTFmut)],2,as.numeric)
  
  library(reshape); library(ggplot2) ;library(scales);library(RCurl);  library (grid)
    
  ### ceRNA target expression
  data.cTarExp = data[,c('sample', ftr.cTarExp)]
  data.cTarExp$variable = rep(ftr.cTarExp, NROW(data)); colnames(data.cTarExp) = c("sample", 'value', 'variable')
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

  ## ceRNA target CNV
  data.cTarCNV = data[,c('sample', ftr.cTarCNV)]
  data.cTarCNV = melt(data.cTarCNV, id.vars='sample')
  data.cTarCNV$value <- factor(sign(as.numeric(as.character(data.cTarCNV$value))), levels = c(-1,0,1))
  data.cTarCNV <- orderSampleLable(data.cTarCNV)
  head(data.cTarCNV)
  
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
    p <- ggplot(data.cTarExp, aes(variable, sample) ) + geom_tile(aes(fill = rescale), color = 'white') 
    p1 <- p + scale_y_discrete(breaks=NULL) + 
              scale_fill_gradient2(low='#99CC00',mid='white',high='#FF6600') + 
              theme(panel.background = element_blank(),
                    axis.title=element_blank(), 
                    axis.text.y = element_blank(),
      #               axis.text.x = element_blank(),            
                    legend.position="none",
      #               axis.title.x = element_text(face="bold", colour="#990000", size=16), axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
                    axis.title.x = element_text(face="bold", colour="#990000", size=14), axis.text.x  = element_text(angle=90, vjust=0.5, size=14))


#     p <- ggplot(data.cRegExp.m, aes(variable, sample)) +   geom_tile(aes(fill = value,height = .2),color = 'white')
#     p2 <- p + scale_y_discrete(breaks=NULL) + 
#               scale_fill_gradient2(low='green',mid='white',high='orange') + 
#               theme_grey(base_size = base_size) + 
#               theme(panel.background = element_blank(),
#                     axis.title=element_blank(), 
#                     axis.text.y = element_blank(),
#       #             axis.text.x = element_blank(), 
#                     legend.position="none" ,
#       #             axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
#                     axis.title.x = element_text(face="bold", colour="#990000", size=14), axis.text.x  = element_text(angle=90, vjust=0.5, size=14))

    p <- ggplot(data.cTarCNV, aes(variable, sample)) +   geom_tile(aes(fill = value ), color = 'white')
    p2 <- p +  scale_y_discrete(breaks=NULL) +  
              scale_fill_manual(values=c("#99CCFF","white","#FF99CC")) +  
              theme(panel.background = element_blank(),
                    axis.title=element_blank(), 
                    axis.text.y = element_blank(),
                    #axis.text.x = element_blank(), 
                    legend.position="none" ,
                    #axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
                    axis.title.x = element_text(face="bold", colour="#990000", size=14), axis.text.x  = element_text(angle=90, vjust=0.5, size= 14))


    p <- ggplot(data.cRegCNV, aes(variable, sample)) +   geom_tile(aes(fill = value ), color = 'white')
    p3 <- p + scale_y_discrete(breaks=NULL) +  
            scale_fill_manual(values=c("#99CCFF","white","#FF99CC")) +  
#           theme_grey(base_size = base_size) +   
            theme(panel.background = element_blank(),
                  axis.title=element_blank(), 
                  axis.text.y = element_blank(),
                    legend.position="none", 
#                     axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#333333", size=10), axis.text.x  = element_text(angle=90, vjust=0.5, size=10))
      
    p <- ggplot(data.cRegTFmut.m, aes(variable, sample)) +  geom_tile(aes(fill = value), colour =   "white") 
    p4 <- p + scale_fill_manual(values=c('white',"#CC3366")) +
            scale_y_discrete(breaks=NULL) + 
            theme(panel.background = element_blank(),
                  axis.title=element_blank(), 
                  axis.text.y = element_blank(), 
                    legend.position="none",
#             axis.title.x = element_text(face="bold", colour="#990000", size=16),axis.text.x  = element_text(angle=90, vjust=0.5, size=16),
            axis.title.x = element_text(face="bold", colour="#990000", size=8), axis.text.x  = element_text(angle=90, vjust=0.5, size=8))
    
     
  len1 = length(unique(data.cTarExp$variable)); len2 = length(unique(data.cRegExp.m$variable)); len3= length(unique(data.cRegCNV$variable)); len4=length(unique(data.cRegTFmut.m$variable)); 
  lensum = len3 + len4
  layt<-grid.layout(nrow=1,ncol=4,widths=c(1/20,len2/lensum * 19/20,len3/lensum * 19/20,len4/lensum * 19/20))

  layt<-grid.layout(nrow=1,ncol=4,widths=c(1/22,1/22, len3/lensum * 10/22,len4/lensum * 10/22))


  pdf(paste(figd,"/heatmap_nose_regCerna_Aug24_",tgene, ".pdf",sep=""),width=15,height=18)
  #   grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(p1,  vp = viewport(layout.pos.row=1,layout.pos.col=1))
  print(p2, vp = viewport(layout.pos.row=1,layout.pos.col=2))
  print(p3, vp = viewport(layout.pos.row=1,layout.pos.col=3))
  print(p4, vp = viewport(layout.pos.row=1,layout.pos.col=4))

  dev.off()

}

calKS = function(inputfile){
  tgene = tail(unlist(strsplit(inputfile,"_")),1)
  if( ! file.exists(inputfile) | file.info(inputfile)$size == 0 ){
    return(rep(NA,5))
  }
  print(inputfile)
  data = read.csv(inputfile, stringsAsFactors= F, header=T,sep="\t")
#   colnames(data ) = c("btname","tgene", "mod1.mse", "mod1.r2",
#                       "mod2.mse", "mod2.r2",
#                       "mod3.mse", "mod3.r2",
#                       "rss1", "rss2", "rss3",
#                       "cntTarCNV", "cntRegCNV", "cntRegTFMut", "cntSmp")
#   
  ks.12 = ks.test(as.numeric(data$mod1.mse), as.numeric(data$mod2.mse),alternative='less')
  ks.23 =  ks.test(as.numeric(data$mod2.mse), as.numeric(data$mod3.mse),alternative='less')
  ks.13 =  ks.test(as.numeric(data$mod1.mse), as.numeric(data$mod3.mse),alternative='less')
  out= sprintf("%s\t%5.3f\t%5.3e\t%5.3f\t%5.3e\t%5.3f\t%5.3e" ,tgene, ks.12$statistic, ks.12$p.value, 
               ks.23$statistic, ks.23$p.value, ks.13$statistic, ks.13$p.value)
  return(out)
}

#### main------
require(scales)
inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/runAug22/" , sep="")
output  = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/summary_runAug22", sep="")

inputfileLst = paste(inputDir, "/", list.files(inputDir,pattern="result"), sep="")
# tgene = 'ABCA1';
# inputfile = paste(inputDir, '/', 'result_', tgene, sep="")
resDF = sapply(inputfileLst, calKS)
resDF = t(resDF)

header <- "tgene\tks12.mse\tks12.mse.pval\tks23.mse\tks23.mse.pval\tks13.mse\tks13.mse.pval"
write(header, file=output)
write(resDF, file=output, append=T)

resDF = read.table(output, sep="\t",stringsAsFactors=F,header=T)
pcut = 0.01

NROW(resDF[which(as.numeric(resDF$ks12.mse.pval) <= pcut),])/ NROW(resDF)
NROW(resDF[which(as.numeric(resDF$ks23.mse.pval) <= pcut),])/ NROW(resDF)
NROW(resDF[which(as.numeric(resDF$ks13.mse.pval) <= pcut),])/ NROW(resDF)

NROW(resDF[which(as.numeric(resDF$ks12.mse.pval) <= pcut),])
NROW(resDF[which(as.numeric(resDF$ks23.mse.pval) <= pcut),])
NROW(resDF[which(as.numeric(resDF$ks13.mse.pval) <= pcut),])
NROW(resDF)


resDF.sig = resDF[which(  as.numeric(resDF$ks12.mse.pval) <= pcut & 
                          as.numeric(resDF$ks23.mse.pval) <= pcut & 
                          as.numeric(resDF$ks13.mse.pval) <= pcut),]


resDF.sig = resDF.sig[order(resDF.sig$ks23.mse.pval),]

head(resDF.sig)
###----plot density
inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/runAug22" , sep="")
plotMSEDensity = function(tgene){
  #   tgene = tail(unlist(strsplit(inputfile,"_")),1)
  inputfile = paste(inputDir, '/', 'result_', tgene, sep="")
  
  if( ! file.exists(inputfile) | file.info(inputfile)$size == 0 ){
    return(rep(NA,5))
  }
  print(inputfile)
  data = read.csv(inputfile, stringsAsFactors= F, header=T,sep="\t")
#   colnames(data ) = c("btname","tgene", "mod1.mse", "mod1.r2",
#                       "mod2.mse", "mod2.r2",
#                       "mod3.mse", "mod3.r2",
#                       "rss1", "rss2", "rss3",
#                       "cntTarCNV", "cntRegCNV", "cntRegTFMut", "cntSmp")
  data[,-c(1,2)] <- apply(data[,-c(1,2)], 2, as.numeric)
  ymax = max(density(as.numeric(data$mod1.mse))$y, density(as.numeric(data$mod2.mse))$y, density(as.numeric(data$mod3.mse))$y)
  xmax = max(data$mod1.mse, data$mod2.mse,data$mod3.mse)
  xmin = min(data$mod1.mse, data$mod2.mse, data$mod3.mse)
  pdf(paste(figd,"/density_",tgene, "_Aug22.pdf", sep=""))
  

    plot(0, type="n", ylab = "density", xlab="",xlim = c(xmin,xmax),
       ylim=c(0,ymax),
       main = paste(tgene, "MSE"))

  lines(density(as.numeric(data$mod1.mse)), col= '#FF9900', lwd = 2)
  lines(density(as.numeric(data$mod2.mse)), col = '#006633', lwd=2) 
  lines(density(as.numeric(data$mod3.mse)), col = '#990033', lwd=2) 
  
  legend('topright',legend=c("model1", "model2", "model3"),lty=1, 
                col=c("#FF9900","#006633","#990033"), bty ='n', cex = .75, lwd=2)
  
    dev.off()
}

plotMulti_MseDensity = function(tglist){
#   pdf(paste(figd,"/density_multiple_Aug24.pdf", sep=""))
  par(mfrow=c(4,4),mar=c(0.1,.1,.1,.1))
  for (tgene in tglist){
    #   tgene = tail(unlist(strsplit(inputfile,"_")),1)
    inputfile = paste(inputDir, '/', 'result_', tgene, sep="")
    
    #     if( ! file.exists(inputfile) | file.info(inputfile)$size == 0 ){
    #       return(rep(NA,5))
    #     }
    print(inputfile)
    data = read.table(inputfile, stringsAsFactors= F, header=T, sep="\t")
    summary(data)
    data[,-c(1,2)] <- apply(data[,-c(1,2)], 2, function(x){as.numeric(as.character(x))})
    
    ymax = max(density(as.numeric(data$mod1.mse))$y, density(as.numeric(data$mod2.mse))$y, density(as.numeric(data$mod3.mse))$y)
    xmax = max(data$mod1.mse, data$mod2.mse,data$mod3.mse)
    xmin = min(data$mod1.mse, data$mod2.mse, data$mod3.mse)
    
    plot(0, type="n", xlim = c(xmin,xmax), xaxt = 'n', yaxt = 'n',
         ylim=c(0,ymax),tcl = 0  )
    axis(side=1,labels=F,tick=F); axis(side=2, labels=F,tick=F)
    lines(density(as.numeric(data$mod1.mse)), col= '#FF9900', lwd= 2)
    lines(density(as.numeric(data$mod2.mse)), col = '#006633', lwd=2) 
    lines(density(as.numeric(data$mod3.mse)), col = '#990033', lwd=2) 
    
    legend('topright',legend=c(tgene, "model 1", "model 2", "model 3"),lty=1, 
           col=c("white", "#FF9900","#006633","#990033"), bty ='n', cex = .75, lwd=2)
  }
  
#   dev.off()
}

plotMSEDensity("SV2A")

inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/runAug22" , sep="")
op = par()
tglist = sample(resDF.sig$tgene,min(NROW(resDF.sig),16))
tglist = as.character(resDF.sig$tgene[1:16])
pdf(paste(figd,"/density_multiple_Aug24_sig.pdf", sep=""))
plotMulti_MseDensity(tglist)
dev.off()

inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/data/" , sep="")
for (tgene in tglist ){
  noseHeatmap(tgene)
}



resDF.sig = resDF[which(  as.numeric(resDF$ks12.mse.pval) > pcut | 
                            as.numeric(resDF$ks23.mse.pval) > pcut | 
                            as.numeric(resDF$ks13.mse.pval) > pcut),]

resDF.nosig = resDF.nosig[order(resDF.nosig$ks13.mse.pval,decreasing=T),]

tglist = resDF.nosig$tgene[1:16]
tglist = as.character(tglist)

inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/result/runAug22" , sep="")
pdf(paste(figd,"/density_multiple_Aug24_nosig.pdf", sep=""))
plotMulti_MseDensity(tglist)
dev.off()

inputDir = paste(rootd, "/DATA/projFocus/result/07152014/reg/tfCerna/data/" , sep="")
for (tgene in tglist ){
  noseHeatmap(tgene)
}



###### 
resDF.sig.plot <- resDF.sig
x = as.numeric(as.character(resDF.sig.plot$ks12.mse.pval)); x = ifelse(x ==0, 1e-200, x)
y = as.numeric(as.character(resDF.sig.plot$ks23.mse.pval)); y = ifelse(y ==0, 1e-200, y)
z = as.numeric(as.character(resDF.sig.plot$ks13.mse.pval)); z = ifelse(z ==0, 1e-200, z)

scatter3D(x =-log10(x), y = -log10(y),  z = -log10(z), pch= 19)
