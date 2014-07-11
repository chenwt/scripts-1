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



permSum = function(infile) {

dataDF = read.table(infile, header = T, sep="\t", stringsAsFactors=F)
real = dataDF[1,]
perm = dataDF[-1,]
tgene = real$tgene
###--look at the correlation
# plot(sort(perm$optCorr,decreasing=T))
# points(real$optCorr,col="red", pch=19)
# points(perm$optCorr[perm$optCorr > real$optCorr], col = "blue")
# 

print (length(perm$optCorr[perm$optCorr > real$optCorr]) / 100)

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
names(yvalue) <- dataDF$tgene
xvalue = -log10(permCorr_Pval[,2])
xvalue[which(xvalue=='NA')] <- 0 
names(xvalue) <- dataDF$tgene

### plot significance level plot
plot(xvalue,yvalue,pch=16,col=mycolor, font = 2, font.lab = 2, xlim=c(1,max(xvalue)), ylim=c(1,max(yvalue)),
     xlab="-log10(pval_perm)", ylab="-log10(pval_full)",main= paste(tgene))
abline(h=2,col = "green", lwd=3)
abline(v=2,col= "green", lwd=3)
legend(x=1.5, y=16,legend=c("negative","positive"),fill=c("blue","red"))
# text(x=2.1,y=8,labels=tgene, font= 2)
points(xvalue[1],yvalue[1], pch = 21, col = colors()[54], lwd = 8)


finalSig = xvalue[xvalue>=2 & yvalue>=2]

print(paste(paste("bothsig",length(names(xvalue[xvalue>=2 & yvalue>=2])),sep=":"),
                  paste("bothnosig",length(xvalue[xvalue<2 & yvalue<2]), sep=":"), 
                  paste("permsig",length(xvalue[xvalue>=2 & yvalue<2]),sep=":"), 
                  paste("fullsig",length(xvalue[xvalue<2 & yvalue>=2]),sep=":"), 
                  paste("finalSig(-1,0]:(0,1]", paste(table(cut(dataDF[names(finalSig),11],breaks=c(-1.1,0,1.1))),sep=":"),sep=":"),
                  #               paste("finalSig(-1,0]:(0,1]", paste(table(cut(dataDF[,11],breaks=c(-1,0,1))),":"),":"), 
                  collapse="\n") )

}


infile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul4/flex_max_1000/permSummary/summary_AP3B1_result_flex_max_1000.tsv"
permSum(infile)

infile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul4/flex_max_1000/permSummary/summary_LIN28B_result_flex_max_1000.tsv"
permSum(infile)

infile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul4/flex_max_1000/permSummary/summary_PPFIA2_result_flex_max_1000.tsv"
permSum(infile)

name = c("summary_UBASH3B_result_flex_max_1000.tsv", "summary_PTGER2_result_flex_max_1000.tsv", "summary_PHF13_result_flex_max_1000.tsv", "summary_PAQR3_result_flex_max_1000.tsv", 
         "summary_LIN28B_result_flex_max_1000.tsv", "summary_FAT1_result_flex_max_1000.tsv", "summary_AP3B1_result_flex_max_1000.tsv","summary_ACSL6_result_flex_max_1000.tsv")

pdf(paste(figd,"/step3-4_greedyOptPermuSummary_07101014", ".pdf",sep=""))

for (x in name){
infile = paste("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul4/flex_max_1000/permSummary/",x,sep="")
permSum(infile)
}

dev.off()
