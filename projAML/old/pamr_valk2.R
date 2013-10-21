setwd("Z:/AML/valk2/")
load("valk2_pamr.Rdata")
require(pamr)
install.packages("matrixStats")
install.packages("bioDist")
source("http://www.bioconductor.org/biocLite.R")
biocLite("bioDist")
require(useful)

expDF <- as.data.frame(expData)
names(clustDefined)<-c("cluster","probeset")
clust <- list()
for( i in 1:16){
ls <- as.character(clustDefined[which(clustDefined$cluster==i),2])
lsAnno <- annoData$ProbeSet
clust[[i]] <-unlist(annoData[unlist(lapply(ls,function(x) grep(x,lsAnno))),1])
}

annoDF <- as.data.frame(merge(annoData,clustDefined,by.x="ProbeSet",by.y="probeset"))
undef <- colnames(expDF[,setdiff(names(expDF)[-(1:2)],annoDF$ID_REF)])
undefDF <- as.data.frame(cbind(ID_REF= undef,cluster= rep(0,length(undef))))
sampleDesign <- rbind(annoDF[,2:3],undefDF)
expDF <- expDF[,sort(names(expDF))]
sampleDesign <- sampleDesign[order(sampleDesign$ID_REF),]
colnames(expDF) <- c(names(expDF)[1:2],as.numeric(sampleDesign$cluster))
subset(expDF,names(expDF)!= "0")

exp.data <- expDF[,which((colnames(expDF)) != "0")]
dataPM <- list(x=as.matrix(exp.data[,-(1:2)]),y=as.numeric(colnames(exp.data)[-(1:2)]),genenames=exp.data[,2],geneid=exp.data[,1])
resT <- pamr.train(dataPM)

names(clustDefined)<-c("cluster","probeset")
expDF <- as.data.frame(expData)
expDF[1:3,setdiff(names(expDF),undef)]
exp270<- expDF[,setdiff(names(expDF),undef)]

sampleDesign <- annoDF[,2:3]
sampleDesign <- sampleDesign[order(sampleDesign$ID_REF),]
exp271 <- exp270
rm("exp270")
dataPM <- list(x=as.matrix(exp271[,-(1:2)]),y=sampleDesign$cluster ,genenames=exp271[,2], geneid=exp271[,1])

#### pamr train process
resT <- pamr.train(dataPM)
resCV <- pamr.cv(fit=resT,data=dataPM,nfold=10)
pamr.plotcv(khan.results)
pamr.plotcv(resCV)
resCV
pamr.confusion(resT, threshold=1.392)
pamr.confusion(resT, threshold=2.783)
pamr.confusion(resT, threshold=4.175)
pamr.plotcvprob(resT, dataPM, threshold=2.783)
pamr.plotcen(resT, dataPM, threshold=2.783)
pamr.plotcenpamr.geneplot(resT, dataPM, threshold=8.350 )
pamr.geneplot(resT, dataPM, threshold=8.350 )
par(mfrow=c(1,1),mar=c(0.3,3,3,3))
pamr.geneplot(resT, dataPM, threshold=8.350 )
resFdr<- pamr.fdr(resT,dataPM)
resFdr<- pamr.fdr(resT,dataPM)
resFdr<- pamr.fdr(resT,dataPM)

centroids.output <- cbind(expDF[,1:2],resT$centroids)

