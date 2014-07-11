#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
rm(list=ls())
usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
typeSelect = "max" #which type of output <all> to select a cutoff <max> to select random optimal(n = 100), or a number to set cutoff value
typeTol = "flex" # < fix> to set tol == 0 , < flex> to select tol adaptively
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


args = getArgs()
if(length(args) < 1 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop(paste(error,"wrong input parameter!"))
}
setwd(system("pwd",intern=T))
input     = args[['greedy']]
param     = args['param']
output = paste(resultflexfile,".significant.summary",sep="")


#---test on local--
input = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul8/optCorr.result_flex_max_1000.tsv"
output = paste(input,".significant.summary",sep="")
param = "flex_max_1000"
#--test-end-

dataDF = read.delim(input, sep="\t", header=T)

# dataDF = dataDF[match(unique(dataDF$tgene), dataDF$tgene),]
rownames(dataDF) <- dataDF$tgene


fullCorr_Pval <- t(vapply(dataDF[,9], FUN=function(x){vapply(unlist(strsplit(gsub("< ","",x),":")),as.numeric,0.1)}, c(0.1,0.1)))
rownames(fullCorr_Pval) <- dataDF[,9]

permCorr_Pval <- t(vapply(dataDF[,10], FUN=function(x){vapply(unlist(strsplit(gsub("< ","",x),":")),as.numeric,0.1)}, c(0.1,0.1)))
rownames(permCorr_Pval) <- dataDF[,10]


source(paste(rootd, "/scripts/myR/jingGraphic.R",sep=""))
dataDF$optCorr <- as.numeric(dataDF$optCorr)
mycolor = val2col(z=dataDF$optCorr,zlim=c(min(dataDF$optCorr),max(dataDF$optCorr)), col=colorRampPalette(c("blue","white","red"))(length(dataDF$optCorr)))

yvalue = -log10(fullCorr_Pval[,2])
yvalue[which(yvalue=='NA')] <- 0
names(yvalue) <- dataDF$tgene
xvalue = -log10(permCorr_Pval[,2])
xvalue[which(xvalue=='NA')] <- 0 
names(xvalue) <- dataDF$tgene

### plot significance level plot
pdf(paste(figd,"/step3-4_greedyOptSummary_07101014_", param, ".pdf",sep=""))
plot(xvalue,yvalue,pch=16,col=mycolor, font = 2, font.lab = 2,
     xlab="-log10(pval_perm)", ylab="-log10(pval_full)",main= paste("Greedy algorithm result \n [", param,"]"))
abline(h=2,col = "green", lwd=3)
abline(v=2,col= "green", lwd=3)
legend(x=1.5, y=16,legend=c("negative","positive"),fill=c("blue","red"))
dev.off()


finalSig = xvalue[xvalue>=2 & yvalue>=2]

fileConn = file(paste(output,".counts",sep=""))
writeLines( paste(paste("bothsig",length(names(xvalue[xvalue>=2 & yvalue>=2])),sep=":"),
              paste("bothnosig",length(xvalue[xvalue<2 & yvalue<2]), sep=":"), 
              paste("permsig",length(xvalue[xvalue>=2 & yvalue<2]),sep=":"), 
              paste("fullsig",length(xvalue[xvalue<2 & yvalue>=2]),sep=":"), 
              paste("finalSig(-1,0]:(0,1]", paste(table(cut(dataDF[names(finalSig),11],breaks=c(-1.1,0,1.1))),sep=":"),sep=":"),
#               paste("finalSig(-1,0]:(0,1]", paste(table(cut(dataDF[,11],breaks=c(-1,0,1))),":"),":"), 
              collapse="\n"), fileConn)

close(fileConn)
# text(x=0,y=1,labels="both not sig",pos=2)
# text(x=3,y=1,labels="full not sig",pos=1)

####------second part
datafinalD <- dataDF[names(finalSig),]
datafinalD <- datafinalD[datafinalD[,11]>0, ]
datafinalD <- datafinalD[order(datafinalD$optCorr,decreasing=T),]

write.table(datafinalD, output,row.names=F,sep="\t",quote=F)
regSmpD = sapply(datafinalD$optGene.optSmp, FUN=function(x){
    paste(t(sapply(unlist(strsplit(gsub("\\]","",gsub("\\[","",as.character(x))),";")), FUN=function(y){
    unlist(strsplit(y,":"))
  }))
})

### output netfile
tarRegDF = t(sapply(1:NROW(datafinalD), FUN=function(i){
    x = datafinalD$optGene.optSmp[i]
    tgene = as.character(datafinalD$tgene[i])
    c(tgene, paste(unique(t(sapply(unlist(strsplit(gsub("\\]","",gsub("\\[","",as.character(x))),";")), FUN=function(y){
    unlist(strsplit(y,":"))
  }))[,1]), collapse=";"))
}))

write.table(tarRegDF, paste(output,".netfile",sep=""),row.names=F,sep="\t",quote=F)

allRegulators = NA
for (i in 1:nrow(datafinalD)){
  allRegulators = c( allRegulators, unlist(strsplit(as.character(tarRegDF[i,2]),";")))
}
allRegulators = unique(allRegulators[-1])

write.table(allRegulators, paste(output,".regulators",sep=""),col.names=F, row.names=F,sep="\t",quote=F)



#######----------plot 
## permutation correlation v.s greedy correlation
permCorr =  vapply(datafinalD$permuCorr.pval, FUN=function(x){as.numeric(unlist(strsplit(as.character(x),":"))[1])}, 1.0)
par(mfrow=c(1,1))
pdf(paste(figd,"/step3-4_greedyOptJuly82014.CorrPermVSOpt_07102014.pdf",sep=""))
plot(datafinalD$optCorr,permCorr, xlim=c(0,1),ylim=c(0,1), pch=19, lwd = 0.1,  col = colors()[563], font = 2,
     xlab = "Optimized correlation", ylab = "Random permutation correlation")
abline(a= 0, b= 1, col = "gray")
dev.off()
########-----------exploratory analysis
# rect(xleft=0,xright=2,ybottom=0,ytop=2,border="green",lwd=2)
# points(x=dataDF$optCorr, y=-log10(permCorr_Pval[,2]), col="orange", pch = 16)
# abline(h = 2, col="red",lwd=2)

####----
#regulator
plot(density(datafinalD$mutReg/datafinalD$totalReg), col="red", lwd=3)
lines(density(datafinalD$optCorrReg/datafinalD$mutReg),col="blue",lwd=3)
lines(density(datafinalD$optCorrReg/datafinalD$totalReg),col="green",lwd=3)
legend(x=0.25, y=7, legend=c("mut_candidate / candidate","driver / mut_candidate","driver / candidate"), 
#        box.lwd = 0,box.col = "white",bg = "white",
       bty="n",
       col=c("red", "blue","green"),lty=1, lwd=2)

plot(datafinalD$mutReg/datafinalD$totalReg, datafinalD$optCorrReg/datafinalD$mutReg, 
     pch=16, col="blue")

#intact sample
plot(density(datafinalD$optCorrSmp/datafinalD$totalSmp),col="green",lwd=3, type='l', lty=2)
lines(density(datafinalD$mutSmp/datafinalD$totalSmp), col="red", lwd=3, type='l', lty=2)
lines(density(datafinalD$optCorrSmp/datafinalD$mutSmp),col="blue",lwd=3, type='l', lty=2)

plot(datafinalD$mutSmp/datafinalD$totalSmp, datafinalD$optCorrSmp/datafinalD$mutSmp, 
     pch=19, col="blue")

###---
pdf(paste(figd,"/step3-4_summary_boxplot_05212014.pdf",sep=""))
boxplot(cbind(allSample=datafinalD[,2], mutatedSample=datafinalD[,3] ,selectedSample=datafinalD[,4]), col="lightblue",
        main="sample number")

# boxplot(cbind(mutatedSample=datafinalD[,3]/datafinalD[,2] * 100 ,selectedSample=datafinalD[,4]/datafinalD[,2] * 100),
#         main="Percentage of samples size")

boxplot(cbind(mutatedReg=datafinalD[,6] ,selectedReg=datafinalD[,7]), col="lightblue",
        main="regulator number")
dev.off()


###---
# datafinalD[order(datafinalD[,6],decreasing=T),1:11][1:100,]
boxplot(cbind(mutatedSample=datafinalD[,6]/datafinalD[,5] * 100 ,selectedSample=datafinalD[,7]/datafinalD[,5] * 100),
        main="Percentage of regulator number")

plot(density(datafinalD$optCorrReg))
lines(density(datafinalD$mutReg),col="blue")
