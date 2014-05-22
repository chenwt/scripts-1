#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript

usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
typeSelect = "max" #which type of output <all> to select a cutoff <max> to select random optimal(n = 100), or a number to set cutoff value
typeTol = "fix" # < fix> to set tol == 0 , < flex> to select tol adaptively
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
figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/May2014/fig/"

source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))


resultfixfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_fix_max_1000.tsv"
resultflexfile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runMay5/optCorr.result_flex_max_1000.tsv"
output = paste(resultflexfile,".significant.summary",sep="")
# 
# resultfile <- resultfixfile
# loadData <- function (resultfile) {
#   
#   con <- file(resultfile, open="r")
#   dataDF = as.data.frame(matrix(NA, ncol=13))
#   colnames(dataDF) = paste(rep("col",13), 1:13, sep="")
#   i = 0
#   while(length((linecrt <- readLines(con,n=1,warn = F) )>0)){
#     i = i+ 1
#     vectcrt <- unlist(strsplit(linecrt,"\t"))
#     if(length(vectcrt)!=13)print(c(i,vectcrt))
#     names(vectcrt) <- colnames(dataDF)
#     dataDF <- rbind(dataDF, vectcrt)
#   }
#   close(con)
#   dataDF <- na.omit(dataDF)
# #   colnames(dataDF) <- dataDF[1,]
# #   dataDF <- dataDF[-1,]
# #   
#   return(dataDF)
# }

fixDF = read.delim(resultfixfile, sep="\t", header=T)
flexDF = read.delim(resultflexfile, sep="\t", header=T)

# dataDF <- fixDF; param <- "fix_max_1000"
dataDF <- flexDF; param <- "flex_max_1000"

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

pdf(paste(figd,"/step3-4_greedyOptSummary_05212014", param, ".pdf",sep=""))
plot(xvalue,yvalue,pch=16,col=mycolor, font = 2, font.lab = 2,
     xlab="-log10(pval_perm)", ylab="-log10(pval_full)",main= paste("Greedy algorithm result \n [", param,"]"))
abline(h=2,col = "green", lwd=3)
abline(v=2,col= "green", lwd=3)
legend(x=1.5, y=16,legend=c("negative","positive"),fill=c("blue","red"))
dev.off()

c(bothsig=length(names(xvalue[xvalue>=2 & yvalue>=2])),bothnosig= length(xvalue[xvalue<2 & yvalue<2]), permsig=length(xvalue[xvalue>=2 & yvalue<2]), fullfig=length(xvalue[xvalue<2 & yvalue>=2]) )

finalSig = xvalue[xvalue>=2 & yvalue>=2]
table(cut(dataDF[names(finalSig),11],breaks=c(-1,0,1)))
table(cut(dataDF[,11],breaks=c(-1,0,1)))

# text(x=0,y=1,labels="both not sig",pos=2)
# text(x=3,y=1,labels="full not sig",pos=1)

####------second part
datafinalD <- dataDF[names(finalSig),]
datafinalD <- datafinalD[datafinalD[,11]>0, ]
datafinalD <- datafinalD[order(datafinalD$optCorr,decreasing=T),]

write.table(datafinalD, output,row.names=F,sep="\t",quote=F)
# paste(datafinalD$tgene,collapse=";")

# rect(xleft=0,xright=2,ybottom=0,ytop=2,border="green",lwd=2)
# 
# points(x=dataDF$optCorr, y=-log10(permCorr_Pval[,2]), col="orange", pch = 16)
# abline(h = 2, col="red",lwd=2)
legend(x=1.5, y=16,legend=c("negative","positive"),fill=c("blue","red"))
pdf(paste(figd,"/step3-4_summary_boxplot_05212014.pdf",sep=""))
boxplot(cbind(allSample=datafinalD[,2], mutatedSample=datafinalD[,3] ,selectedSample=datafinalD[,4]), col="lightblue",
        main="sample number")

# boxplot(cbind(mutatedSample=datafinalD[,3]/datafinalD[,2] * 100 ,selectedSample=datafinalD[,4]/datafinalD[,2] * 100),
#         main="Percentage of samples size")

boxplot(cbind(mutatedReg=datafinalD[,6] ,selectedReg=datafinalD[,7]), col="lightblue",
        main="regulator number")
dev.off()

datafinalD[order(datafinalD[,6],decreasing=T),1:11][1:100,]

boxplot(cbind(mutatedSample=datafinalD[,6]/datafinalD[,5] * 100 ,selectedSample=datafinalD[,7]/datafinalD[,5] * 100),
        main="Percentage of regulator number")

plot(density(datafinalD$optCorrReg))
lines(density(datafinalD$mutReg),col="blue")
