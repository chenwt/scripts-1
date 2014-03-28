##---init---
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
figd = paste(rootd,"/DATA/projFocus/report/Mar2014/fig/",sep="")
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))
setwd(jxy(rootd, "/DATA/projFocus/data/03102014/tcgal2som"))


file = "genome.wustl.edu__Illumina_All.maf.matrix.promoter2k.Mar-14-2014.matrix"
fdata = read.delim2(file,skip=1)
# dim(fdata)
fmut = fdata[,-c(1,2,3,4)]
fmut = unique(fmut)[,-c(1,2,3)]
# dim(fmut)
frecur = rowSums(fmut)

pdf(jxy(figd,"/barplot_Mutspectrum_",cdt,".pdf"))
barplot(table(frecur)/length(frecur)*100,col="lightblue",width=1,space=0.2,
        ylab="frequency", xlab = "recurence(N = 1001)",
        main = "brca promoters region mutation recurrence(+/-2k)",
        ylim=c(0,100))
text(x=seq(0.5,7.5,1) + 0.2,y=20,
     labels=round(table(frecur)/length(frecur),3)*100,font=2,cex=0.8)

file = "genome.wustl.edu__Illumina_All.maf.matrix.utr3p.Mar-14-2014.matrix"
fdata = read.delim2(file,skip=1)
fmut = fdata[,-c(1,2,3,4)]
fmut = unique(fmut)[,-c(1,2,3)]
# dim(fmut)
frecur = rowSums(fmut)
barplot(table(frecur)/length(frecur)*100,col="lightblue",width=1,space=0.2,
        main = " brca 3' UTR region mutation recurrence",
        ylab="frequency", xlab = "recurence(N = 1001)", ylim=c(0,100))
text(x=seq(1:length(table(frecur))) - 0.5 + 0.2,y=20,
     labels=round(table(frecur)/length(frecur),3)*100,font=2,cex=0.8)

dev.off()

file = "genome.wustl.edu__Illumina_All.maf.matrix.utr5p.Mar-14-2014.matrix"
fdata = read.delim2(file,skip=1)
fmut = fdata[,-c(1,2,3,4)]
fmut = unique(fmut)[,-c(1,2,3)]
# dim(fmut)
frecur = rowSums(fmut)
table(frecur)


