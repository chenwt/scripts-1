
require(useful)
require("pamr")

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
sampleDesign <- annoDF[,2:3]

expDF <- expDF[,sort(names(expDF))]
sampleDesign <- sampleDesign[order(sampleDesign$ID_REF),]
# colnames(expDF) <- c(names(expDF)[1:2],as.numeric(sampleDesign$cluster))

exp271<- expDF[,setdiff(names(expDF),undef)]

dataPM <- list(x=as.matrix(exp271[,-(1:2)]),y=sampleDesign$cluster ,genenames=exp271[,2], geneid=exp271[,1])
resT <- pamr.train(dataPM)
resCV <- pamr.cv(fit=resT,data=dataPM,nfold=10)
pamr.plotcv(resCV)
pamr.confusion(resT, threshold=1.392)
pamr.plotcvprob(resT, dataPM, threshold=2.783)

pamr.plotcen(resT, dataPM, threshold=2.783)
pamr.geneplot(resT, dataPM, threshold=8.350 )

resFdr<- pamr.fdr(resT,dataPM)
pamr.plotfdr(resFdr)
pamr.listgenes(resT, exp271, threshold=8.35)

