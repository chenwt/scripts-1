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
figd = paste(rootd,"/DATA/projFocus/report/topDown_02042014/fig/",sep="")
cdt = paste(unlist(strsplit(date()," "))[c(2,3,5)],collapse="-")

##----sec1 barPlot for 
#73 WGS samples
gsFreqF = paste(rootd,"/DATA/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_geneSampleNumber.txt",sep="")
gsFreq = read.table(gsFreqF,sep="\t",header=F,stringsAsFactors=T)
sf = gsFreq[order(gsFreq[,2],decreasing=T),2]
st = gsFreq[1,]
sf = as.numeric(sf)
sfCut = table(cut(sf,breaks=c(seq(9,max(sf),by=10),max(sf))))

pdf(paste(figd,"barplot_focusTargets",cdt,".pdf",sep=""))
par(mar=c(6,3,3,0),mgp=c(4.2,0.5,0))
barplot(sfCut/length(sf) * 100,
        las=2,col = "lightblue",
        xlab="Number of samples",
        ylab = "Percentage",
        main="Distribution of Gint sample count"
)
text(x=6,y=20,cex=0.8, font = 2,
     labels=paste("Total # of Targets: ", 321 ,"\nTotal # of Samples: 73",sep="")
)

regFreqF = paste(rootd,"/DATA/projFocus/result/02022014/geneSamples/brca_gslist_combCG_gintact_Mar-7-2014.txt.deg_20140307.txt_regulatorFreq.txt",sep="")
regFreq = read.table(regFreqF,stringsAsFactors=F,header=T,skip=1)
rf = regFreq[order(regFreq[,2],decreasing=T),]
gt = rf$gene
rf = as.numeric(rf[,2])
table(cut(rf,breaks=c(-1,0,1,2,3,4,5,10,730)))
names(rf) = gt

# pdf(paste(figd,"barplot_focusTargets",cdt,".pdf",sep=""))
par(mar=c(6,3,3,0),mgp=c(4.2,0.5,0))
barplot(table(cut(rf,breaks=c(0,seq(10,300,by=50),400, 500,max(rf)))) / length(rf) * 100,
        las=2,col = "lightblue",
        xlab="Number of CeRNET Regulators",#ylab = "Freq of CeRNET Targets"
        main="Distribution of Regulators numbers for Targets",
        )
text(x=8,y=25,cex=0.8, font = 2,
     labels=paste("Total # of Targets: ",length(rf[rf>0]),sep="")
     )
dev.off()

####0----somatic mutation
fsom = "som.mat"
wd          = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/CHEK1-temp/"
dataRaw = read.delim2(fsom)
smps_temp = vapply(colnames(dataRaw),FUN=function(x){substr(x,6,15)},'a')
colnames(dataRaw) = smps_temp  
gMut = dataRaw[,1:4]
gene = 'KIF23'
idx = (grep(gene,gMut$ame,perl=T))
dataSom = apply(subset(dataRaw[idx,],select=smps),2,as.numeric  )
dataSom = dataSom[which(rowSums(dataSom)>0),]

mycol = colorRampPalette(c("white","red"))(25)

figd = paste(rootd,"/DATA/projFocus/report/topDown_02042014/fig/",sep="")

pdf(paste(figd,"/heatmap_CHEK1_KIF23.pdf",sep=""))
par(mar=c(8,0,2,0))
heatmap.2(dataSom, col=mycol,trace='none',dendrogram="none",cexRow=1,
          labRow="\n\nGenomic\nlocation\ngroups\n(1kb)", 
          key=F, 
          sepwidth=c(0.05, 0.05),  # width of the borders
          sepcolor='black', scale="none",
          main=paste("somatic mutation of",gene),
          )

sf2 = paste(rootd,"/DATA/projFocus/result/02022014/som/test/CHEK1_KIF23.som.mat",sep="")

dataRaw = read.delim2(sf2)
smps_temp = vapply(colnames(dataRaw),FUN=function(x){substr(x,6,15)},'a')
colnames(dataRaw) = smps_temp  
gMut = dataRaw[,1:4]
gene = 'KIF23'
idx = (grep(gene,gMut$ame,perl=T))
dataSom = apply(subset(dataRaw[idx,],select=smps),2,as.numeric  )
dataSom = dataSom[which(rowSums(dataSom)>0),]

figd = paste(rootd,"/DATA/projFocus/report/topDown_02042014/fig/",sep="")
pdf(paste(figd,"/heatmap_CHEK1_KIF23_dist0.pdf",sep=""))
par(mar=c(8,0,2,0))
heatmap.2(dataSom, col=mycol,trace='none',dendrogram="none",cexRow=1,
          labRow="\n\nGenomic\nlocation\ngroups\n(point)", 
          key=F, 
          sepwidth=c(0.05, 0.05),  # width of the borders
          sepcolor='black', scale="none",
          main=paste("somatic mutation of",gene),
)
dev.off()


###-------------------mutations freqs
pdf(paste(figd,"hist_mutationFreqPerGene",cdt,".pdf",sep=""))
par(mfrow=c(1,2),mar=c(3,3,3,0),mgp=c(1.5,0.5,0))
ff = paste(rootd,"/DATA/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist0_Mar-7-2014.matrix.MutFreqForGene",sep="")
dfreq = read.table(ff)
hist(dfreq$V2,breaks=c(seq(0,max(dfreq$V2),by=100),max(dfreq$V2 + 1)),freq=T,
     xlab = "number of mutations", 
     main="point mutation numbers per gene")
ff = paste(rootd,"/DATA/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist1000_Mar-7-2014.matrix.MutFreqForGene",sep="")
dfreq = read.table(ff)
hist(dfreq$V2,breaks=c(seq(0,max(dfreq$V2),by=50),max(dfreq$V2 + 1)),freq=T,
     xlab = "number of mutation hotspots",
     main = "mutation hotspots per gene")
dev.off()

pdf(paste(figd,"hist_mutationOcurrenceInSample",cdt,".pdf",sep=""))
par(mfrow=c(1,2),mar=c(3,3,3,0),mgp=c(1.5,0.5,0))
ff = paste(rootd,"/DATA/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist0_Mar-7-2014.matrix.mutFreq",sep="")
dfreq = read.table(ff)
hist(dfreq$V5,
#      breaks=c(seq(0,max(dfreq$V5),by=100),max(dfreq$V5 + 1)),freq=T,
     xlab = "mutations occurence", 
     main="mutation penetriablity")

ff = paste(rootd,"/DATA/projFocus/result/02022014/som/brca_somTumorWU_combinedCG_Regulator_v3_Mar-6-2014.mat_dist1000_Mar-7-2014.matrix.mutFreq",sep="")
dfreq = read.table(ff)
hist(dfreq$V5,
#      breaks=c(seq(0,max(dfreq$V5),by=5),max(dfreq$V5 + 1)),
     xlab = "mutation hotspots occurence",
     main = "mutation hotspots penetriablity ")
dev.off()
