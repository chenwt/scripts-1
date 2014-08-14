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

noseHeatmap = function(inputdata){
  data = inputdata
  data = data[order(as.numeric(as.character(data$exp))),]
  allregs = setdiff(colnames(data),c('exp','cnv'))
  data = data[, c("exp", "cnv", sort(allregs))]
  
  data$sample = row.names(data)
  colnames(data)
  library(reshape); library(ggplot2) ;library(scales);library(RCurl);  library (grid)
  dataTfMut  = data[,c('sample',allregs)]
  dataTfMut.m = melt(dataTfMut)
  dataTfMut.m$value <- factor(dataTfMut.m$value, levels=c(0,1))
  dataTfMut.m <- orderSampleLable(dataTfMut.m)
  
  sum(as.numeric(as.character(dataTfMut[,3])))
  
  dataCNV = data.frame(cbind(sample=data$sample,variable = rep('      cnv',NROW(data)), value=data$cnv))
  dataCNV$value <- factor(sign(as.numeric(as.character(dataCNV$value))), levels = c(-1,0,1))
  dataCNV <- orderSampleLable(dataCNV)
  head(dataCNV)
  
  dataExp = data.frame(cbind(sample=data$sample,variable = rep('      exp',NROW(data)), 
                             value=data$exp))
  dataExp$value <- as.numeric(as.character(dataExp$value))
  dataExp <- orderSampleLable(dataExp)
  dataExp$rescale <- rescale(dataExp$value,to=c(-1,1))
  
  
    p <- ggplot(dataTfMut.m, aes(variable, sample)) +
        geom_tile(aes(fill = value), colour =   "white") 
    pmut <- p + scale_fill_manual(values=c('white',"red")) + 
        theme(axis.title=element_blank(), 
          axis.text.y = element_blank(),
          legend.position="none",
          axis.title.x = element_text(face="bold", colour="#990000", size=10),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=10))
    pmut
  
    p <- ggplot(dataCNV, aes(variable,sample)) + 
        geom_tile(aes(fill = value), color = 'white')
    pcnv <- p + scale_fill_manual(values=c("lightblue","white",colors()[119])) +  
      theme(axis.title=element_blank(), 
          axis.text.y = element_blank(),
          legend.position="none",
          axis.title.x = element_text(face="bold", colour="#990000", size=16),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=16))
  
    p <- ggplot(dataExp, aes(variable,sample)) + 
          geom_tile(aes(fill = rescale), color = 'white')
    pexp <- p + scale_fill_gradient2(low='green',mid='white',high='orange') + 
      theme(axis.title=element_blank(), 
            axis.text.y = element_blank(),
            legend.position="none" ,
            axis.title.x = element_text(face="bold", colour="#990000", size=16),
            axis.text.x  = element_text(angle=90, vjust=0.5, size=16))
  
  
  
  pdf(paste(figd,"/heatmap_diagnois_1_Aug12_",tgene, ".pdf",sep=""))
  layt<-grid.layout(nrow=1,ncol=3,widths=c(1/8,1/8,6/8))
  #   grid.show.layout(layt)
  grid.newpage()
  pushViewport(viewport(layout=layt))
  print(pexp,  vp = viewport(layout.pos.row=1,layout.pos.col=1))
  print(pcnv, vp = viewport(layout.pos.row=1,layout.pos.col=2))
  print(pmut, vp = viewport(layout.pos.row=1,layout.pos.col=3))
  dev.off()
}

tarTFmutReg = function(tgene, plot = "False"){
  inputdata= paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/",tgene,".tgTFreg.input",sep="")
  inputdata2 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_",tgene,".temp",sep="")
  inputdata3 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_act_",tgene,"_1000_1000_mutSampleVector",sep="")
  
  mutSmp = read.table(inputdata3, header=T,stringsAsFactors=F); smpMerg = colnames(mutSmp)[which(mutSmp[1,]>0)]
  mutD = read.table(inputdata2,header=T,stringsAsFactors=F); dataExp = t(mutD[1,]); mutD = data.frame(t(mutD[-1,] )); 
  allsmp1 = rownames(mutD); mutD[setdiff(allsmp1,smpMerg),] <- 0
  length(which(rowSums(mutD)>0))
  
  datacnv = t(read.table(inputdata, header=T,stringsAsFactors=F))
  head(datacnv)
  allsmp2 = rownames(datacnv)
  allsmp = intersect(allsmp1, allsmp2); datacnv <- datacnv[allsmp,]
  
  data = data.frame(cbind(dataExp, mutD)); tgene = colnames(data)[1]; allregs = colnames(data)[-1]; colnames(data) <- c('exp',unlist(as.character(allregs)))
  data$cnv  <- 0; 
  data[allsmp,'cnv'] <-  unlist(as.numeric(datacnv[,2]))
  
  data[setdiff(allsmp,smpMerg),allregs] <-0
  
  ftest = anova(lm(exp~cnv,data),lm(exp~.,data[,c('exp','cnv',allregs)]))
  res = c(tgene = tgene,  formatTestOut(ftest),allregs = paste(allregs, collapse=","))
  if (plot == 'True'){noseHeatmap(data)}
  return(res)
}

###----main-------
gfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/summary/runAug13_0.05_fail.genelist"
tgenelist = unlist(read.table(gfile,stringsAsFactors=F))


tgene = tgenelist[1]
tgene = 'OMG'
require(scales)
testNormality = function(tgene){
    ### exp-cnv--mut(wrong)
    inputdata= paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/",tgene,".tgTFreg.input",sep="")
    ### exp-mut(beforeGreedy)
    inputdata2 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_",tgene,".temp",sep="")
    ### target-collapsedTFMut(afterGreedy)
    inputdata3 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/data/tfGreedy_act_",tgene,"_1000_1000_mutSampleVector",sep="")
    print(tgene)
    if (!(file.exists(inputdata) && file.exists(inputdata2) && file.exists(inputdata3))) {
      return(c('NA','NA','NA'))
    }
    mutSmp = read.table(inputdata3, header=T,stringsAsFactors=F); smpMerg = colnames(mutSmp)[which(mutSmp[1,]>0)]
    mutD = read.table(inputdata2,header=T,stringsAsFactors=F); dataExp = t(mutD[1,]); mutD = data.frame(t(mutD[-1,] )); 
    allsmp1 = rownames(mutD); mutD[setdiff(allsmp1,smpMerg),] <- 0
    length(which(rowSums(mutD)>0))
    dataExp = rescale(dataExp,to=c(1,10))
    res1 = c(tgene = tgene, before=shapiro.test(dataExp)$p.value, log2=shapiro.test(log2(dataExp))$p.value)
#     return(res1)
    datacnv = t(read.table(inputdata, header=T,stringsAsFactors=F))
    head(datacnv)
    allsmp2 = rownames(datacnv)
    allsmp = intersect(allsmp1, allsmp2); datacnv <- datacnv[allsmp,]
    
    data = data.frame(cbind(dataExp, mutD)); tgene = colnames(data)[1]; allregs = colnames(data)[-1]; colnames(data) <- c('exp',unlist(as.character(allregs)))
    data$cnv  <- 0; 
    data[allsmp,'cnv'] <-  unlist(as.numeric(datacnv[,2]))
    
    data[setdiff(allsmp,smpMerg),allregs] <-0
    
    data$exp <- log2(data$exp)
    ftest = anova(lm(exp~cnv,data),lm(exp~.,data[,c('exp','cnv',allregs)]))
    res2 = c(formatTestOut(ftest),allregs = paste(allregs, collapse=","))
    
    return(c(res1,res2))
}

resultLog2 = vapply(tgenelist, FUN=testNormality,rep('a',10))

resultLog2 = t(resultLog2)
write.table(file=paste(gfile,".log2_result",sep=""), x= resultLog2, col.names=T, row.names=F, quote = F, sep="\t")
head(resultLog2)


plot(log10(as.numeric(resultLog2[3,])))

targenes = unlist(read.table("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug13/input.genelist.data",stringsAsFactors=F))
result <- rep(NA,8 )
for (tgene in targenes){
  tgene
  result <- rbind(result,tarTFmutReg(tgene))
}

result = result[-1,]; result <- result[order(result[,7]),]
write.table(result,file="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/summary/reg_tfTarMut_runAug13.txt",
            col.names=T,row.names=T, quote=F, sep="\t")

inputdata= paste("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/",tgene,".tgTFreg.input",sep="")
data = read.table(inputdata, header=T,stringsAsFactors=F)

plot(density(log2(rescale(as.numeric(data[1,]),to=c(1,5)))))
plot(density(rescale(as.numeric(data[1,]),to=c(1,5))))
