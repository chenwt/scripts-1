###Rscript
##J.HE
##Generate plots for report on Marww2 lab meeting

setRootd = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysnameß']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/projFocus/"
  }
  return(rootd)
}

rootd = setRootd()
CDT = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")

# numReg = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethFree.10smapMore.deg_20140325.txt.10more_stat.GintRegCount"
# numSmp = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_Mar-24-2014_CnvMethFree.10smapMore.deg_20140325.txt.10more_stat.GintSmpCount"
# output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslist_Mar-24-2014_CnvMethFree_10more_Gint_stats_"
# totalSmp = 778
# 

numSmp = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_gint_deg_Mar-25-2014_stat.GintSmpCount"
numReg = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/gslist/gslist_gint_deg_Mar-25-2014_stat.GintRegCount"
output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Mar2014/fig/gslist_Mar-24-2014_CnvMethSomFree_10more_Gint_stats_"
totalSmp = 778

##----sec1 barPlot for 
library(xtable)
require("Hmisc")

gsFreq = read.table(numSmp,sep="\t",header=F,stringsAsFactors=T)
sf = gsFreq[order(gsFreq[,2],decreasing=T),2]
st = gsFreq[1,]
sf = as.numeric(sf)
sfCut = table(cut(sf,breaks=c(seq(9,max(sf),by=50),max(sf))))
sfCut.table = xtable(as.data.frame(round(sfCut/length(sf) * 100, digits=2)))
dvips(latex(sfCut.table),file=paste(output, "_table_smpCount_", CDT, ".ps",sep=""))

regFreq = read.table(numReg,stringsAsFactors=F,header=T,skip=1)
rf = regFreq[order(regFreq[,2],decreasing=T),]
gt = rf$gene
rf = as.numeric(rf[,2])
names(rf) = gt
rfCut = table(cut(rf,breaks=c(0,seq(10,300,by=50),400, 500,max(rf))))
rfCut.table = xtable(as.data.frame(round(rfCut/length(rf) * 100,digits=2)))
dvips(latex(rfCut.table),
      file=paste(output, "_table_regCount_", CDT, ".ps",sep=""))

pdf(paste(output,"_barplot", CDT, ".pdf", sep=""))
par(mar=c(6,3,3,0),mgp=c(4.2,0.5,0))
barplot(sfCut/length(sf) * 100,
        las=2,col = "lightblue",
        xlab="Number of samples",
        ylab = "Percentage",
        main="Distribution of Gint sample count"
)
text(x=2,y=20,cex=0.8, font = 2,
     labels=paste("Total # of Targets: ", length(rf[rf>0]) ,"\nTotal # of Samples: ", totalSmp ,sep="")
)



par(mar=c(6,3,3,0),mgp=c(4.2,0.5,0))
barplot( rfCut / length(rf) * 100,
        las=2,col = "lightblue",
        xlab="Number of CeRNET Regulators",#ylab = "Freq of CeRNET Targets"
        main="Distribution of Regulators numbers for Gint targets",
        )
text(x=8,y=25,cex=0.8, font = 2,
     labels=paste("Total # of Targets: ",length(rf[rf>0]),sep="")
     )
dev.off()

