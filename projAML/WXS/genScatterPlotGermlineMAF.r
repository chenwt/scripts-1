##Desp: used in projAML, WXS, germline mutation number, maf visulization
#input: input data(raw counts of specified germline coordinates)
#output: scatterplot 1: relapse v.s. tumore mutaions number;
#        scatterplot 2: relaspe v.s. tumor maf( 4 * 4 layout)
sysInfo = Sys.info()
ifelse(sysInfo['sysname'] == "Darwin",
       setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result_anno_filter_final/germlineFinal/"),
       setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/callVars/result_anno_filter_final/germlineFinal/")
)


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
dataRef =  apply(data,c(1,2),function(x){
  temp = as.character(x)
  ref = as.numeric(unlist(strsplit(temp,"/"))[2])
  return(ref)
})
rm(data)

colnames(dataAlt) = patientID
dataAltN = dataAlt[,seq(3,50,by=3)]
dataAltT = dataAlt[,seq(2,49,by=3)]
dataAltR = dataAlt[,seq(1,48,by=3)]
colnames(dataAltN) =patientID[seq(3,50,by=3)]
colnames(dataAltT) =patientID[seq(2,49,by=3)]
colnames(dataAltR) =patientID[seq(1,48,by=3)]
rm(dataAlt)

dataRefN = dataRef[,seq(3,50,by=3)]
dataRefT = dataRef[,seq(2,49,by=3)]
dataRefR = dataRef[,seq(1,48,by=3)]
colnames(dataRefN) =patientID[seq(3,50,by=3)]
colnames(dataRefT) =patientID[seq(2,49,by=3)]
colnames(dataRefR) =patientID[seq(1,48,by=3)]
rm(dataRef)

##---raw data processing
mafN = dataAltN / (dataAltN  + dataRefN)
mafT = dataAltT / (dataAltT  +  dataRefT)
mafR = dataAltR / (dataAltR  + dataRefR)
patientID = sapply(colnames(mafN),function(x){substr(x,0,6)})
save(mafR,mafT,mafN,patientID,file="germlineMAF.rda")

require(RColorBrewer)
pdf("scatterPlotGermlineMAF.pdf")
layout(matrix(1:16,byrow=T,nrow=4))
par(mar=c(2.8,2.8,0.5,0.5),mgp=c(1.5,0.5,0))
for (indexPt in 1:16){
  color = rep("gray",length(mafT[,indexPt]))
  x = mafN[,indexPt]; y = mafT[,indexPt]; z = mafR[,indexPt]
  indexZero = intersect(which(x ==0),intersect(which(y == 0), which(z==0)))
  x = x[-indexZero];y = y[-indexZero];z = z[-indexZero]
  color = rep("blue",length(x))
  xylim = c(0,1)
  plot(x=z,y=y,
     col=color,
     xlab = paste(patientID[indexPt],"Re"), ylab=paste(patientID[indexPt],"Tu"),
     xlim =xylim, ylim=xylim,
     cex=0.2)
}
dev.off()

#------cluster before plotting
#patient PANVGP 
indexPt = 8
no = mafN[,indexPt];tu = mafT[,indexPt]; re = mafR[,indexPt]
indexZero = intersect(which(no ==0),intersect(which(tu == 0), which(re==0)))
no = no[-indexZero]; tu = tu[-indexZero]; re = re[-indexZero]
dmat = as.matrix(cbind(no,tu,re))
dmat = na.omit(dmat)

## hirarchical clustering...
ddist = dist(dmat)
save(mafN,mafR,mafT,ddist,file="germlineMAFDist.rda")
hfit = hclust(ddist)
pdf("plotHClustGermline_v1.pdf")
plot(fit)
dev.off()


## svm
require(e1071)
fit = svm(dmat)

##kmeans-- not working 
fit = kmeans(dmat,centers=3)
plot(dmat[,-1],col=fit$cluster,cex=0.3)

dcls2 = dmat[which(fit$cluster==1,arr.ind=T),]
fit2 = kmeans(dcls2,centers=3)
plot(dcls2[,-1],col=fit2$cluster,cex=0.3)

##Gaussian Mixture
install.packages("mclust")
require(mclust)
mclust = Mclust(dmat,G=3)
save(mclust,file="clustGM.rda")
##density estimation DBSCAN
install.packages("fpc")
require()

