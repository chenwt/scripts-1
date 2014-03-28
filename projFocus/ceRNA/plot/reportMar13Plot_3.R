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
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

wd          = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/CHEK1-temp/"
figd = paste(rootd,"/DATA/projFocus/report/topDown_02042014/fig/",sep="")
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")

require(gplots)
mf1 = paste(rootd,"/DATA/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist0_Mar-7-2014.matrix.mutFreq.FreqNum",sep="")
mutFreq1 = read.table(mf1)

mutFreq1t = table(mutFreq1)
pdf(paste(figd,"barplot_mutationRecurrence_",cdt,".pdf",sep=""),height=10)
par(mar=c(6,4,0,0))
barplot(/nrow(mutFreq1) * 100,width=1.5,space=0.2,
        ylim = c(0,100),
        col="lightblue",border="gray",las=2,
        xlab = "mutation hotspot recurrence",
        ylab = "percentage")

dev.off()


mf1k = paste(rootd, "/DATA/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist1000_Mar-7-2014.matrix.mutFreq",sep="")
mutFreq1k = read.table(mf1k)[,5]

pdf(paste(figd,"barplot_mutationHotspotRecurrence_",cdt,".pdf",sep=""),height=10)
par(mar=c(6,4,0,0))
barplot(table(mutFreq1k[-1])/length(mutFreq1k) * 100,width=1.5,space=0.2,
     ylim = c(0,20), 
     col="lightblue",border="gray",las=2,
     xlab = "mutation hotspot recurrence",
     ylab = "percentage")
axis()
dev.off()

