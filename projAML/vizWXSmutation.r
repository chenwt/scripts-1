#!/bin/Rscirpts
#Author: Jing He
#Date: Jun 6
#Last Update: Jun 6
#parameters: input files: /ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/annotated
#example: Rscript 

###------------------------------header------------------------------
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters")
	exit
}else{
	print(args)
}

###------------------------------start coding##------------------------------
install.packages("reshape2")
setwd("/Volumes//ac_lab/jh3283/SCRATCH/projAML/WXS/annotated/")
fns <- dir()

dataRaw <- data.frame()
for (i in 1:length(fns)){
  temp <- read.table(fns[i],header=F,stringsAsFactors=F,sep="\t")[,-15]
  dataRaw <- rbind(dataRaw,
                cbind(PID=rep(gsub(".filtered.eff.vcf.parsed...var","",fns[[i]]),times=nrow(temp)),
                      temp))
}

data <- dataRaw[,c(1:3,9,11:12)]
colnames(data) <- c("PID","Chr","Pos","Gene","snpEff","Eff")

#filterring

dataCoding <- data[which(data$snpEff=="START_GAINED"|data$snpEff=="STOP_GAINED"|data$snpEff=="NON_SYNONYMOUS_CODING"),]
require(reshape2)
dataCod.cast <- dcast(data=dataCoding,PID~Gene,value.var="PID",fun.aggregate=length)
dataGene <- as.matrix(dataCod.cast[,-1])
dataGene[dataGene>1] <- 1
row.names(dataGene) <- dataCod.cast[,1]
table(dataGene)

op <- par(mar = c(5,5,4,2) + 0.1)
image(dataGene,col=c("white","red"),axes=F)
axis(side = 2, labels = rownames(dataGene), 
     at = seq(0, by = 0.5, length.out = nrow(dataGene)), las=1)
axis(side = 1, labels = colnames(dataGene),at=seq(0,by=0.5,length.out=ncol(dataGene)),las=1)
box()
par(op)

install.packages("bipartite")
library(bipartite)
mat<-dataGene
rownames(mat)<-rownames(dataGene)
colnames(mat)<-colnames(dataGene)
visweb(mat,type="None",labsize=2,square="b",box.col="red") 

plot(y=dataCoding$Chr,x=dataCoding$Pos,pch="*",ylab = "Pos", xlab="Chr")
for (i in 1:11)  abline(h=i,col="blue")