##Desp: used in projAML, WXS, germline mutation number, maf visulization
#input: input data(raw counts of specified germline coordinates)
#output: scatterplot 1: relapse v.s. tumore mutaions number;
#        scatterplot 2: relaspe v.s. tumor maf( 4 * 4 layout)
sysInfo = Sys.info()
ifelse(sysInfo['sysname'] == "Darwin",
       setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result_anno_filter_final/germlineFinal/"),
       setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result_anno_filter_final/germlineFinal/")
)

##-----mutation number
data = read.table("finalTSV/germlineMutationCount.txt",header=T,sep="\t")
corTest = cor.test(data$Relapse,data$Tumor)
paste("correlation:",round(corTest$estimate,3)," pvalue:",round(corTest$p.value,3),sep="")

X = data$Relapse
Y = data$Tumor
label= data$PID
pdf("scatterPlotGermlineCount.pdf")
par(mar=c(2.8,2.8,2.8,2.8),mgp=c(1.5,0.5,0))
plot(x=data$Relapse,y=data$Tumor,
     pch=15,col="blue",
     xlim = c(min(min(data$Tumor),min(data$Relapse)),max(max(data$Tumor),max(data$Relapse))), ylim = c(min(min(data$Tumor),min(data$Relapse)),max(max(data$Tumor),max(data$Relapse))),
     xlab="Relapse", ylab="Tumor")
text(x=X,y = Y - 0.01 * max(Y),labels=label,
     cex = 0.5)
abline(a=0,b=1,col="gray")
dev.off()

##-----for MAF
data = read.table("germlineMutationAlleleFreqRawcountSorted.txt",header=F,sep="\t")
data = t(data)
patientID = data[1,][-1]
chrPos = data[,1][-1]
data = data[,-1]
data = data[-1,]

dataAlt =  apply(data,c(1,2),function(x){
    temp = as.character(x)
    alt = as.numeric(unlist(strsplit(temp,"/"))[1])
    return(alt)
})
colnames(dataAlt) = patientID
dataAltN = dataAlt[,seq(3,50,by=3)]
dataAltT = dataAlt[,seq(2,49,by=3)]
dataAltR = dataAlt[,seq(1,48,by=3)]

dataRef =  apply(data,c(1,2),function(x){
  temp = as.character(x)
  ref = as.numeric(unlist(strsplit(temp,"/"))[2])
  return(ref)
})
colnames(dataRef) = patientID
dataRefN = dataRef[,seq(3,50,by=3)]
dataRefT = dataRef[,seq(2,49,by=3)]
dataRefR = dataRef[,seq(1,48,by=3)]

##---raw data processing
indes = intersect(which(dataAltT[,1] > 1),which(dataRefT[,1] > 1))

mafN = dataAltN[indes,] / dataRefN[indes,]
mafT = dataAltT[indes,] / dataRefT[indes,]
mafR = dataAltR[indes,] / dataRefR[indes,]

indes = intersect(which(dataAltT[,1] > 1),which(dataRefT[,1] > 1))
  color = rep("gray",length(mafT[,1]))
  y = mafT[,1]
  x = mafR[,1]
  index = intersect(which(x>= 0.8),which(y <= 0.2 ))
  color[index] = "green"
  index = intersect(which(x<= 0.1),which(y >= 0.9 ))
  color[index] = "yellow"

  
mycolor = color()
#xylim = c(0,max(mafR[indes,1],mafT[indes,1] ))
xylim = c(0,1)
plot(x=mafR[,1],
     y=mafT[,1],
     xlab="Relapse",ylab="Tumor",
#     xlim = xylim, ylim = xylim,     
     col=color,
     cex=0.3)
points(log10(mafN[,1]),log10(mafR[,1]),col="blue")

colNames = sapply(patientID,function(x){substr(x,0,6)})
###plot MAF as x-y plot for each sample
germlineMAFbySampleType = function(){
  layout(matrix(1:16,byrow=T,nrow=4))
  par(mar=c(2.8,2.8,0.5,0.5),mgp=c(1.5,0.5,0))
  for (i in seq(1,48,by=3)){
  print(colNames[(i+2)])
  plot(x=log10(dataAlt[,(i+2)]), y=log10(dataRef[,(i+2)]),
       pch="*",col=rgb(0,0,0,alpha=0.5),cex=0.5,
       xlab=paste(colNames[(i+2)],"Alt"), ylab= paste(colNames[(i+2)],"Ref"),
    )
  points(x=log10(dataAlt[,(i+1)]),y=log10(dataRef[,(i+1)]),pch="*",col=rgb(0,1,0,alpha=0.5),cex=0.5)
  points(x=log10(dataAlt[,i]),y=log10(dataRef[,i]),pch="*",col=rgb(1,0,0,alpha=0.5),cex=0.5)  
  }
}

pdf("scatterPlotFinalMAFGermlineAll.pdf")
germlineMAFbySampleType()
dev.off()


#pdf("scatterPlotFinalMAFGermlineReplase.pdf")
#germlineMAFbySampleType(type=0)
#dev.off()


rawcountPlot = apply(data,c(1,2),function(x){
  temp = as.character(x)
  alt = as.numeric(unlist(strsplit(temp,"/"))[1])
  ref = as.numeric(unlist(strsplit(temp,"/"))[2])
  return(alt/ref)
})

###generate color matrix
#for (i in seq(1,48,by=3)){
  if(rawcountPlot[])
}
  

pdf("scatterPlotFinalMAFGermline.pdf")
layout(matrix(1:16,byrow=T,nrow=4))
par(mar=c(2.8,2.8,0.5,0.5),mgp=c(1.5,0.5,0))
for (i in seq(1,48,by=3)){
  print(paste(i+2,i+1))
  plot(rawcountPlot[c(i,i+1),] * 100, 
       pch="*",col="blue",cex=1.5,
       xlim=c(0,100),ylim=c(0,100)
  )
}
dev.off()


