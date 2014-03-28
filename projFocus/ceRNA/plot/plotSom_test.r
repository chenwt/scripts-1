
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


###---set---
cwd = paste(rootd,"DATA/projFocus/result/02022014/som/test",sep="")
setwd(cwd)
require(gplots)
require(useful)


##---data
f = "brca_somTumorWU_APC_Regulator__dist1000.matrix"
dataSom = read.table(f,header=T,sep="\t")
genes  = as.character(unique(dataSom$geneName))
barcodes = vapply(colnames(dataSom)[-c(1:4)],FUN=function(x){substr(x,6,15)},'a')
ds   =   as.matrix(dataSom[,-c(1:4)])
colnames(ds) = barcodes


cumsum(table(cut(sort(table(dataSom$geneName)),breaks=c(0,100,300,500,700,800) )) / length(genes))
par(mar=c(8,4,1,0))
barplot(table(cut(sort(table(dataSom$geneName)),breaks=c(0,100,300,500,700,800) )),
        las=2,cex.axis=0.5,
        sub="number of mutation groups(1kb)", ylab = "number of genes", )
mutGrps = as.character(paste(dataSom$mutChr,dataSom$mutRegStart,dataSom$mutRegEnd))


f = "APC.smps"
dataSmps = unlist(strsplit(readLines(f),"\t|;",perl=T))
target = dataSmps[1]
dataSmps = dataSmps[-1]
dsSub = subset(ds,select=dataSmps)

##---plot
mycol = colorRampPalette(c("white","red"))(25)
pdf(paste(rootd,"/DATA/projFocus/report/plotSom_test_Mar72014.pdf",sep=""))
par(mar=c(10,1,1,0))
heatmap.2(dsSub[grep(genes[1],dataSom$geneName),],col=mycol,trace="none", 
          cexRow=0.3,cexCol=0.5,dendrogram='none'  )
dev.off()
