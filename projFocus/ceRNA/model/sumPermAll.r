##initiate
rm(list=ls())
usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
timeStart = Sys.time()

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
figd = paste(rootd, "/DATA/projFocus/report/Jul2014/fig/",sep="")

source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

##--init end 

##--get cmd argv
args = getArgs()
if(length(args) < 1 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop(paste(error,"wrong input parameter!"))
}
setwd(system("pwd",intern=T))
permTsvDir = args[['dir']]
param     = args['param']
output = paste(resultflexfile,".significant.summary",sep="")

##--test--
permTsvDir = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/sigReg/runGATA3/flex_max_1000/permSummary"
out = "step3-4_greedyOptPermuSummary_GATA3_aracne_07222014"
output = paste(permTsvDir, "/summary_45_07222014.txt", sep="")

##-test-end

##--func---
permSum = function(infile) {
  dataDF = read.table(infile, header = T, sep="\t", stringsAsFactors=F)
  for (i in c(2:8,11)){
    dataDF[,i] = as.numeric(as.character(dataDF[,i]))
  }
  
  tgene = dataDF$tgene[1]
  
  rownames(dataDF) = c(paste(tgene, "actual", sep = "_"), paste(tgene, 1: (NROW(dataDF) -1 ), sep="_"))
  real = dataDF[1,]
  perm = dataDF[-1,]
  

  
  fullCorr = vapply(dataDF$fullCorr.pval,FUN=function(x){as.numeric(unlist(strsplit(x,":"))[1])}, 0.1)
  names(fullCorr) = rownames(dataDF)
  if(fullCorr[1] > 0){
    cntfullCorrBig = length(fullCorr[fullCorr >= fullCorr[1]]) -1 
  }else if (fullCorr[1] < 0 ){
    cntfullCorrBig = length(fullCorr[fullCorr <= fullCorr[1]]) -1 
  }
  
  if(dataDF$totalCorr[1] > 0){
    cntBigtotal = length(fullCorr[fullCorr >= dataDF$totalCorr[1] ]) 
  }else if (dataDF$totalCorr[1] < 0 ){
    cntBigtotal = length(fullCorr[fullCorr <= dataDF$totalCorr[1] ]) 
  }
  
  
  
  if(real$optCorr > 0){
    cntCorrBig = length(perm$optCorr[perm$optCorr >= real$optCorr]) - 1
  }else if (real$optCorr < 0 ){
    cntCorrBig = length(perm$optCorr[perm$optCorr <= real$optCorr]) -1 
  }
  
  
  ###--plot correlations
  par(mfrow=c(1,2), oma = c(0, 0, 3, 0))
  plot(x=1:length(fullCorr), y =fullCorr,pch=19, col= colors()[325],  ylim = c(-1,1), 
       xlab="Randomization", ylab = " All mutation expression correlation ", 
       main = "random mutation correlation ")
  points(x=1, y=fullCorr[1],pch=19,col="red")
  points(as.integer(length(fullCorr)/2), dataDF$totalCorr[1],col = colors()[614], pch =19, lwd = 5)
  
  plot(x = 1: length(dataDF$optCorr), y =dataDF$optCorr,pch=19,  col= colors()[325], ylim = c(-1,1), 
       xlab="Randomization", ylab = " Optimized correlation ", 
       main = "optimized correlation")
  points(x=1, y=dataDF$optCorr[1],pch=19,col="red")
  points(as.integer(length(fullCorr)/2), dataDF$totalCorr[1],col = colors()[614], pch =19, lwd = 5)
  mtext(paste(tgene,":", NROW(dataDF)-1, " correlations", ""), outer = TRUE, cex = 1.5, font =2)
  
  par(mfrow=c(1,1), oma=c(0,0,0,0))

  ## plot -end 
  
  

  ###---look at the pvalue
  
  fullCorr_Pval <- t(vapply(dataDF[,9], FUN=function(x){vapply(unlist(strsplit(gsub("< ","",x),":")),as.numeric,0.1)}, c(0.1,0.1)))
  rownames(fullCorr_Pval) <- dataDF[,9]
  
  permCorr_Pval <- t(vapply(dataDF[,10], FUN=function(x){vapply(unlist(strsplit(gsub("< ","",x),":")),as.numeric,0.1)}, c(0.1,0.1)))
  rownames(permCorr_Pval) <- dataDF[,10]
  
  source(paste(rootd, "/scripts/myR/jingGraphic.R",sep=""))
  dataDF$optCorr <- as.numeric(dataDF$optCorr)
  mycolor = val2col(z=dataDF$optCorr,zlim=c(min(dataDF$optCorr),max(dataDF$optCorr)), col=colorRampPalette(c(colors()[563],"white",colors()[526]))(length(dataDF$optCorr)))
  
  yvalue = -log10(fullCorr_Pval[,2])
  yvalue[which(yvalue=='NA')] <- 0
  names(yvalue) <- rownames(dataDF)
  xvalue = -log10(permCorr_Pval[,2])
  xvalue[which(xvalue=='NA')] <- 0 
  names(xvalue) <- rownames(dataDF)
  
  ### plot significance level plot
  plot(xvalue,yvalue,pch=16,col=mycolor, font = 2, font.lab = 2, xlim=c(1,max(xvalue)), ylim=c(1,max(yvalue)),
       xlab="-log10(pval_perm)", ylab="-log10(pval_full)",main= paste(tgene, ":", NROW(perm)))
  abline(h=2,col = "green", lwd=3)
  abline(v=2,col= "green", lwd=3)
  legend(x=1.5, y=max(yvalue) - 1.5,legend=c("negative","positive"),fill=c("blue","red"))
  # text(x=2.1,y=8,labels=tgene, font= 2)
  points(xvalue[1],yvalue[1], pch = 21, col = colors()[654], lwd = 8)
  
  
  finalSig = xvalue[xvalue>=2 & yvalue>=2]
  cntSigCorrBig = length(which(dataDF[names(finalSig),11] > real$optCorr) == TRUE)
#   print(paste(paste("bothsig",length(names(xvalue[xvalue>=2 & yvalue>=2])),sep=":"),
#                     paste("bothnosig",length(xvalue[xvalue<2 & yvalue<2]), sep=":"), 
#                     paste("permsig",length(xvalue[xvalue>=2 & yvalue<2]),sep=":"), 
#                     paste("fullsig",length(xvalue[xvalue<2 & yvalue>=2]),sep=":"), 
#                     paste("finalSig(-1,0]:(0,1]", paste(table(cut(dataDF[names(finalSig),11],breaks=c(-1.1,0,1.1))),sep=":"),sep=":"),
#                     #               paste("finalSig(-1,0]:(0,1]", paste(table(cut(dataDF[,11],breaks=c(-1,0,1))),":"),":"), 
#                     collapse="\n") )


  
  return(c(tgene = tgene,
           totalperm = as.character(NROW(perm)),
           cntfullCorrBig = as.character(cntfullCorrBig),
           cntBigtotal = as.character(cntBigtotal),
           cntCorrBig = as.character(cntCorrBig), 
           cntSigCorrBig = as.character(cntSigCorrBig), 
           totalCorr = as.character(dataDF$totalCorr[1]),
           actulCorr = as.character(fullCorr[1]),
           actualOptCorr = as.character(dataDF$optCorr[1]),
            bothsig = as.character(length(names(xvalue[xvalue>=2 & yvalue>=2])) ) ,
            bothnosig = as.character( length(xvalue[xvalue<2 & yvalue<2]) ), 
               permsig = as.character( length(xvalue[xvalue>=2 & yvalue<2]) ),
               fullsig = as.character( length(xvalue[xvalue<2 & yvalue>=2]) ),
              finalsig = table(cut(dataDF[names(finalSig),11],breaks=c(-1.1,0,1.1)))
          ))
                
}
##--funcEnd

name = sort(list.files(path=permTsvDir,pattern="*tsv"))
summaryStat = rep(NA,15)
pdf(paste(figd,"/", out, ".pdf",sep=""))
for (x in name){
print(x)
infile = paste(permTsvDir,"/",x,sep="")
summaryStat = rbind(summaryStat, permSum(infile))
}
dev.off()
summaryStat = data.frame(summaryStat[-1,])
rownames(summaryStat) = summaryStat$tgene
summaryStat = summaryStat[order(as.numeric(as.character(summaryStat[,3])),decreasing=F),]




cbind(as.character(summaryStat[,1]), vapply(summaryStat[,3],FUN=function(x){as.numeric(as.character(x))},1) / vapply(summaryStat[,2], FUN=function(x){as.numeric(as.character(x))}, 1) )
cbind(as.character(summaryStat[,1]),summaryStat$actulCorr,summaryStat$totalCorr, 
      vapply(summaryStat$actulCorr,FUN=function(x){abs(as.numeric(as.character(x)))},1) - vapply(summaryStat$totalCorr, FUN=function(x){abs(as.numeric(as.character(x)))}, 1) )

# summaryStat[order(as.numeric(as.character(summaryStat[,8])),decreasing=T),]

write.table(file=output, summaryStat, quote=F,row.names=T,sep="\t")


###-----additional analysis
gata3DriGene = unlist(read.table("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/sigReg/runGATA3/target_drivers_hasGATA3",stringsAsFactors=F))
sumGata3 <- summaryStat[gata3DriGene,]; sumGata3[order(as.numeric(as.character(sumGata3[,3])),decreasing=F),]
cbind(sumGata3[,1], sumGata3[,3], sumGata3[,2], vapply(sumGata3[,3],as.numeric,1) / vapply(sumGata3[,2], as.numeric, 1) )

sumGata3[order(vapply(sumGata3[,3],as.numeric,1) / vapply(sumGata3[,2], as.numeric, 1),decreasing=F),]

sumGata3 = rep(NA,14); cnt = 0 
pdf(paste(figd,"/", out, "_16gata3Drived.pdf",sep=""))
for (x in name){
  if (!is.na(match(unlist(strsplit(x,"_"))[2], gata3DriGene))) {
    cnt = cnt + 1
    infile = paste(permTsvDir,"/",x,sep="")
    sumGata3 = rbind(sumGata3, permSum(infile))
  }
}
data.frame(sumGata3)
sumGata3$
sumGata3 = sumGata3[-1,]; sumGata3[order(as.numeric(as.character(sumGata3[,3])),decreasing=F),]
print(cnt)
dev.off()
