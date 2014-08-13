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

##func---###

formatTestOut = function(test_aov){
  Res.DF_cnv_all = paste(test_aov$Res.Df, collapse="/")
  RSS_cnv_all = paste(round(test_aov$RSS,digits=2),collapse="/")  
  result = c(Res.DF_cnv_all= Res.DF_cnv_all, RSS_cnv_all= RSS_cnv_all, round(unlist(data.frame(test_aov)[2,3:6]),5))
  return(result)
}

anova_cnvVstfmut = function(inputdata){
  
  ##first version
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
      fitlm_all = lm(exp ~ . , data)
      tmp = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
      if(as.numeric(tmp[6])<p_prev){
        res = tmp; p_rev = as.numeric(tmp[6])
      }
    }
  }else{
    
    fitlm_selfcnv = lm(exp ~ cnv, data)
    fitlm_all = lm(exp ~ . , data)
    
    res = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
  }
  return(c(tgene=tgene,cntCNVsmp, cntTfmutSmp,res))
}


outputCoeff = function(inputdata){
  tgene =   unlist(strsplit( tail(unlist(strsplit(x= inputdata, "/")),1), "\\."))[1]
  data = data.frame(t(read.table(inputdata, header=T, sep="\t",row.names=1,stringsAsFactors=F)))
  data = data.frame(apply(data, c(1,2), as.numeric))
  fitlm_selfcnv = lm(exp ~ cnv, data)
  fitlm_all = lm(exp ~ . , data)
  allCoeff = c(fitlm_selfcnv$coefficients[-1], 
               fitlm_all$coefficients[-1])
  mutSmp = which(rowSums(data[,-c(1,2)])>0,arr.ind=T); 
  
  calKS(data$exp[mutSmp],data$exp[-mutSmp])
  res = formatTestOut(anova(fitlm_selfcnv, fitlm_all))
  
  return(c(tgene=tgene,cntCNVsmp, cntTfmutSmp,res))
  
}

orderSampleLable = function(dataExp){
  dataExp$sample <- as.character(dataExp$sample)
  #Then turn it back into an ordered factor
  dataExp$sample <- factor(dataExp$sample, levels=unique(dataExp$sample))
  return(dataExp)
}

noseHeatmap = function(inputdata){
  data = data.frame(inputdata)
#   tgene =   unlist(strsplit( tail(unlist(strsplit(x= inputdata, "/")),1), "\\."))[1]
#   data = data.frame(t(read.table(inputdata, header=T, sep="\t",row.names=1,stringsAsFactors=T)))
  data = data[order(as.numeric(as.character(data$exp))),]
  allregs = setdiff(colnames(data),c('exp','cnv'))
  data = data[, c("exp", "cnv", sort(allregs))]
  
  data$sample = row.names(data)
  colnames(data)
  library(reshape); library(ggplot2) ;library(scales);library(RCurl);  library (grid)
  dataTfMut  = data[,c('sample',allregs)]
  dataTfMut.m = melt(dataTfMut,id.vars='sample')
  dataTfMut.m$value <- factor(dataTfMut.m$value, levels=c(0,1))
  dataTfMut.m <- orderSampleLable(dataTfMut.m)
  
  head(dataTfMut.m)
  
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

### main
# inputDir= "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/"
# allInputFile =  vapply(list.files(path=inputDir,pattern="input"), FUN=function(x){
#   paste(inputDir,"/",x,sep="")
# },'aa')
# 
# result = sapply(allInputFile,anova_cnvVstfmut)
# 
# result = data.frame(t(result));rownames(result) = result[,1]; 
# result = result[,-1]

tarTFmutReg = function(tgene){
  #       tgene = 'FUT8'
  inputdata= paste("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/",tgene,".tgTFreg.input",sep="")
  inputdata2 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/runJuly27/data/tfGreedy_",tgene,".temp",sep="")
  inputdata3 = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/runJuly27/data/tfGreedy_act_",tgene,"_1000_1000_mutSampleVector",sep="")
  
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
  head(data)
  data[setdiff(allsmp,mutSmp),allregs] <-0
  
  ftest = anova(lm(exp~cnv,data),lm(exp~.,data[,c('exp','cnv',allregs)]))
  res = c(tgene = tgene,  formatTestOut(ftest),allregs = paste(allregs, collapse=","))
  print(res)
  noseHeatmap(data)
  return(res)
}

targenes = c("HS3ST1", "IGDCC3", "LBR", "NCALD", "PAPD7", "PPARGC1B", "RAB12", "SLC45A1", "SLC16A6", "AMFR", "USP6NL", "B4GALT5", "CAMLG", "CDKN1B", "ERBB4", "FMNL2", "FSTL4", "FUT8", ## GATA3
             "ATP7A","BMPR1A","NASP","ARL5A")
result <- rep(NA,8 )
for (tgene in targenes){
  
  result <- rbind(result,tarTFmutReg(tgene))
}
###---main----
# inputDir= "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/"
# allInputFile =  vapply(list.files(path=inputDir,pattern="input"), FUN=function(x){
#   paste(inputDir,"/",x,sep="")
# },'aa')
# 
# result = sapply(allInputFile,anova_cnvVstfmut)
# 
# result = data.frame(t(result));rownames(result) = result[,1]; 
# result = result[,-1]
# 
# write.table(result,file="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/summary/reg_tfTarMut_runAug9.txt",
#             col.names=T,row.names=T, quote=F, sep="\t")
# length(which(as.numeric(as.character(result[,8]))<0.05))
# result[(which(as.numeric(as.character(result[,8]))>0.05)),]


####function


###

### test


## LEO1
notar = NA
for (tgene in targenes){
  inputdata= paste("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/reg/tfTar/runAug9/",tgene,".tgTFreg.input",sep="")
  if(file.exists(inputdata)) {
    noseHeatmap(inputdata)
  }else{ 
    notar  = c(notar, tgene)
  }
}
notar <- notar[-1]

###tables of each coefficients
## data


df <- data.frame(a = seq(0, 90, 10), b = seq(10, 100, 10))
df.plot <- ggplot(data = df, aes(x = seq(1, 100, 10))) + 
  geom_line(aes(y = a), colour = 'red') +
  geom_line(aes(y = b), colour = 'blue') +
  scale_x_continuous(breaks = seq(0,100,10))

# make dummy labels for the table content
df$lab <- month.abb[ceiling((df$a+1)/10)]

df.table <- ggplot(df, aes(x = a, y = 0,
                           label = lab, colour = b)) +
  geom_text(size = 3.5) + 
  theme_minimal() + 
  scale_y_continuous(breaks=NULL)+
  theme(panel.grid.major = element_blank(), legend.position = "none",
        panel.border = element_blank(), axis.text.x =  element_blank(),
        axis.ticks =  element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) 

# silly business to align the two plot panels    
gA <- ggplotGrob(df.plot)
gB <- ggplotGrob(df.table)

maxWidth = grid::unit.pmax(gA$widths[2:3], gB$widths[2:3])
gA$widths[2:3] <- as.list(maxWidth)
gB$widths[2:3] <- as.list(maxWidth)

require(gridExtra)
grid.arrange(gA, gB, ncol=1, heights=c(10,1))
install.packages("maps")
install.packages("OIdata");
library(ggplot2)
library(gridExtra)

# line breaks between words for levels of birds$effect:
data(birds)
library(ggplot2)

levels(birds$effect) <- gsub(" ", "\n", levels(birds$effect))
ggplot(birds,
       aes(x = effect,
           y = speed)) +
  geom_boxplot() +   coord_flip()


#### sample size
pdf(paste(figd,"/boxplot_cnvsmpVStfmutsmp_Aug09.pdf",sep=""))
boxplot(as.numeric(as.character(result$V2)), as.numeric(as.character(result$V3)), names=c("cnv_sample","tfmut_sample"), ylab = "event sample size", font=2,pch=19)
dev.off()

NROW(result[(which(as.numeric(as.character(result[,8]))<0.05)),])

kruskal.test(exp ~ cnv, data=data)
kruskal.test(exp~., data=data)
