
CWD ="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/fucFilt/step5-2_dnaseSite/"
setwd(CWD)

###---scan encode 
cgfilelist = list.files(path=CWD,pattern="^cg_brca_greedy_ceRNAdirverMutation.all_06182014.bed.wg")
allfilelist = list.files(path=CWD,pattern="^brca_greedy_ceRNAdirverMutation.all_06182014.bed.wg")

###---cancer genes----
dataDF = as.data.frame(matrix(NA, ncol=7))
for( i in 1:length(cgfilelist)) { 
  infile = paste(CWD,"/",cgfilelist[i],sep="")
  fgroup = gsub(".Jun-18-2014","",
                gsub(pattern="cg_brca_greedy_ceRNAdirverMutation.all_06182014.bed.wgEncodeUwDnaseMcf7",
                     "",cgfilelist[i]))
  fgroup = gsub("Rep1.broadPeak","", fgroup)
  fgroup = gsub("Rep2.broadPeak","", fgroup)
  fgroup = gsub("Rep1.narrowPeak","", fgroup)
  fgroup = gsub("Rep2.narrowPeak","", fgroup)
  
  con <- file(infile, open="r")
  while (length(linecrt <- readLines(con, n=1, warn = F)) > 0){
    vectcrt = c(unlist(strsplit(linecrt, "\t"))[c(4,1,2,3,11,12)],fgroup)
    dataDF <- rbind(dataDF,vectcrt)
  }
  close(con)
}

dataDF = na.omit(dataDF)
colnames(dataDF) = c("gene","chr","mutps","mutpe","score","signalValue","sample")
rownames(dataDF) = paste(dataDF$chr,"-",dataDF$mutps,"-",dataDF$mutpe,paste="")
head(dataDF)

dataDF$signalValue <- as.numeric(dataDF$signalValue)
table(dataDF$sample)
data100HotDF <- dataDF[which(dataDF$sample == "Est100nm1hHotspots"),]
datactrlHotDF <- dataDF[which(dataDF$sample == "Estctrl0hHotspots"),]
datanoneHotDF <- dataDF[which(dataDF$sample == "Hotspots"),]

length(intersect(datanoneHotDF$gene, intersect(datactrlHotDF$gene, data100HotDF$gene)))
require(reshape)

datanoneHotDF2 = cast(datanoneHotDF,formula=gene+chr+mutps+mutpe~sample,value="signalValue",fun.aggregate=mean)
rownames(datanoneHotDF2) = paste(datanoneHotDF2$chr,"-",datanoneHotDF2$mutps,"-",datanoneHotDF2$mutpe,paste="")

datactrlHotDF2  = cast(datactrlHotDF,formula=gene+chr+mutps+mutpe~sample,value="signalValue",fun.aggregate=mean)
rownames(datactrlHotDF2) = paste(datactrlHotDF2$chr,"-",datactrlHotDF2$mutps,"-",datactrlHotDF2$mutpe,paste="")

data100HotDF2   = cast(data100HotDF,formula=gene+chr+mutps+mutpe~sample,value="signalValue",fun.aggregate=mean)
rownames(data100HotDF2) = paste(data100HotDF2$chr,"-",data100HotDF2$mutps,"-",data100HotDF2$mutpe,paste="")

idNewRow = setdiff(rownames(data100HotDF2), c(rownames(datactrlHotDF2), rownames(datanoneHotDF2)))
###--resultNew ---
resultNewSignal <- data100HotDF2[idNewRow,]
resultNewSignal = resultNewSignal[resultNewSignal$Est100nm1hHotspots > 5,]

###--resultDiffExpressed---
idCommCtrlRow = intersect(rownames(data100HotDF2),rownames(datactrlHotDF2))
idCommNoneRow = intersect(rownames(data100HotDF2),rownames(datanoneHotDF2))
t2c.fc  = (data100HotDF2[idCommCtrlRow,5] - datactrlHotDF2[idCommCtrlRow,5]) / datactrlHotDF2[idCommCtrlRow,5]
names(t2c.fc) <- rownames(data100HotDF2[idCommCtrlRow,])
cbind(data100HotDF2[names(t2c.fc[t2c.fc>2]),], t2cFC=t2c.fc[t2c.fc>2])

t2n.fc  = (data100HotDF2[idCommNoneRow,5] - datanoneHotDF2[idCommNoneRow,5]) / datanoneHotDF2[idCommNoneRow,5]
names(t2n.fc) <- rownames(data100HotDF2[idCommNoneRow,])
data100HotDF2[names(t2n.fc[t2n.fc>2]),]

idResult <- setdiff(names(t2c.fc[t2c.fc>2]), names(t2n.fc[t2n.fc>2]))

resultDiffFC <- cbind(data100HotDF2[idResult,], t2cFC=t2c.fc[idResult], t2nFC = t2n.fc[idResult])

###---result final and output--

result = rbind(resultDiffFC, 
               cbind(resultNewSignal, t2cFC=rep(NA, nrow(resultNewSignal)), 
                      t2nFC=rep(NA, nrow(resultNewSignal))))
outfile = paste(CWD,"/mutInDnas/cg_brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv",sep="")
write.table(result,
            file=outfile,
            sep="\t",row.names=F,col.names=T,quote=F
            )

###---all genes-----

###---cancer genes----
filelist = allfilelist
dataDF = as.data.frame(matrix(NA, ncol=7))
for( i in 1:length(filelist)) { 
  infile = paste(CWD,"/",filelist[i],sep="")
  fgroup = gsub(".Jun-18-2014","",
                gsub(pattern="brca_greedy_ceRNAdirverMutation.all_06182014.bed.wgEncodeUwDnaseMcf7",
                     "",filelist[i]))
  fgroup = gsub("Rep1.broadPeak","", fgroup)
  fgroup = gsub("Rep2.broadPeak","", fgroup)
  fgroup = gsub("Rep1.narrowPeak","", fgroup)
  fgroup = gsub("Rep2.narrowPeak","", fgroup)
  
  con <- file(infile, open="r")
  while (length(linecrt <- readLines(con, n=1, warn = F)) > 0){
    vectcrt = c(unlist(strsplit(linecrt, "\t"))[c(4,1,2,3,11,12)],fgroup)
    dataDF <- rbind(dataDF,vectcrt)
  }
  close(con)
}

dataDF = na.omit(dataDF)
colnames(dataDF) = c("gene","chr","mutps","mutpe","score","signalValue","sample")
dataDF$signalValue <- as.numeric(dataDF$signalValue)

table(dataDF$sample)
data100HotDF <- dataDF[which(dataDF$sample == "Est100nm1hHotspots"),]
datactrlHotDF <- dataDF[which(dataDF$sample == "Estctrl0hHotspots"),]
datanoneHotDF <- dataDF[which(dataDF$sample == "Hotspots"),]

require(reshape)

datanoneHotDF2 = cast(datanoneHotDF,formula=gene+chr+mutps+mutpe~sample,value="signalValue",fun.aggregate=mean)
rownames(datanoneHotDF2) = paste(datanoneHotDF2$chr,"-",datanoneHotDF2$mutps,"-",datanoneHotDF2$mutpe,paste="")

datactrlHotDF2  = cast(datactrlHotDF,formula=gene+chr+mutps+mutpe~sample,value="signalValue",fun.aggregate=mean)
rownames(datactrlHotDF2) = paste(datactrlHotDF2$chr,"-",datactrlHotDF2$mutps,"-",datactrlHotDF2$mutpe,paste="")

data100HotDF2   = cast(data100HotDF,formula=gene+chr+mutps+mutpe~sample,value="signalValue",fun.aggregate=mean)
rownames(data100HotDF2) = paste(data100HotDF2$chr,"-",data100HotDF2$mutps,"-",data100HotDF2$mutpe,paste="")

idNewRow = setdiff(rownames(data100HotDF2), c(rownames(datactrlHotDF2), rownames(datanoneHotDF2)))
###--resultNew ---
resultNewSignal <- data100HotDF2[idNewRow,]
resultNewSignal = resultNewSignal[resultNewSignal$Est100nm1hHotspots > 5,]

###--resultDiffExpressed---
idCommCtrlRow = intersect(rownames(data100HotDF2),rownames(datactrlHotDF2))
idCommNoneRow = intersect(rownames(data100HotDF2),rownames(datanoneHotDF2))
t2c.fc  = (data100HotDF2[idCommCtrlRow,5] - datactrlHotDF2[idCommCtrlRow,5]) / datactrlHotDF2[idCommCtrlRow,5]
names(t2c.fc) <- rownames(data100HotDF2[idCommCtrlRow,])
cbind(data100HotDF2[names(t2c.fc[t2c.fc>2]),], t2cFC=t2c.fc[t2c.fc>2])

t2n.fc  = (data100HotDF2[idCommNoneRow,5] - datanoneHotDF2[idCommNoneRow,5]) / datanoneHotDF2[idCommNoneRow,5]
names(t2n.fc) <- rownames(data100HotDF2[idCommNoneRow,])
data100HotDF2[names(t2n.fc[t2n.fc>2]),]

idResult <- setdiff(names(t2c.fc[t2c.fc>2]), names(t2n.fc[t2n.fc>2]))

resultDiffFC <- cbind(data100HotDF2[idResult,], t2cFC=t2c.fc[idResult], t2nFC = t2n.fc[idResult])

###---result final and output--

result = rbind(resultDiffFC, 
               cbind(resultNewSignal, t2cFC=rep(NA, nrow(resultNewSignal)), 
                     t2nFC=rep(NA, nrow(resultNewSignal))))
outfile = paste(CWD,"/mutInDnas/brca_greedy_ceRNAdirverMutation.all.wgEncodeUwDnaseMcf7Est100nm1hHotspots.06182014.tsv",sep="")
write.table(result,
            file=outfile,
            sep="\t",row.names=F,col.names=T,quote=F
)

