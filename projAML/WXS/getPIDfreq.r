#!/usr/bin/Rscript
#Author: J.He
#Date:Sep.15 
#TODO: not finished at all, get error msg when dcast, segfault from C stack overflow

pid <- "PASFEW"
fns <- c(paste(pid,"-NoA.freq",sep=""),
	 paste(pid,"-TuA.freq",sep=""),
	 paste(pid,"-ReA.freq",sep=""))
readData <- function(x) {
  tmp <- read.table(x,header=T,sep="\t")
  tmp <- data.frame(chrPos=paste(tmp$Chr,tmp$Pos,sep="_"),
		    Total=tmp$Totalreads, Alt=tmp$Altreads,
		    pidSample=rep(fns[i],nrow(tmp)) )
  return(tmp)
}
data1 <- readData(fns[1])
data2 <- readData(fns[2])
data3 <- readData(fns[3])

for(i in 1:3) {
  tmp <- readData(fns[i])
  ifelse(i!=1,data <- rbind(data,tmp),data <- tmp)
}

require(reshape2)
nn <- length(unique(data$chrPos))
dataNew <- data.frame(chrPos=unique(data$chrPos),NTotal=rep(0,nn),
		      NAlt=rep(0,nn),TTotal=rep(0,nn),TAlt=rep(0,nn),
		      RTotal=rep(0,nn),RAlt=rep(0,nn))
row.names(dataNew) <- dataNew$chrPos
dataNew[data$chrPos[grep("NoA",data$pidSample)],c(2,3)] <- data[grep("NoA",data$pidSample),c("Total","Alt")]
dataNew[data$chrPos[grep("TuA",data$pidSample)],c(4,5)] <- data[grep("TuA",data$pidSample),c("Total","Alt")]
dataNew[data$chrPos[grep("ReA",data$pidSample)],c(6,7)] <- data[grep("ReA",data$pidSample),c("Total","Alt")]

sapply(unique(dataNew$pidSample),
data.melt <- melt(data,id.vars=c(1,4))
dataNew <- dcast(data.melt,chrPos ~ pidSample + variable, value.var="value")


