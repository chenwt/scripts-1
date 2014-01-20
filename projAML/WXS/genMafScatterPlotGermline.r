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

rawcountPlot = apply(data,c(1,2),function(x){
  temp = as.character(x)
  alt = as.numeric(unlist(strsplit(temp,"/"))[1])
  ref = as.numeric(unlist(strsplit(temp,"/"))[2])
  return(alt/ref)
})

###generate color matrix
for (i in seq(1,48,by=3)){
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


