#!/bin/Rscirpts
#Author: Jing He
#Date: Mar.26,2013
#Last Update: Mar.28,2013
#Modelling the reads using beta distribution

##------------------------------null distribution
source("~/Dropbox/scripts/projNET/testCNV/CNVTest.r")

# n.file <- gsub("CHR","10","/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chrCHR.freq")
snp.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt"

##------------------------------get alpha & beta for null distribution
getParaNull <- function(freq.adj){
	af <- freq.adj$Altfreqadj
	#based on chromosome level
	m <- mean(af)
	v <- var(af)
	if(v < m * (1-m))
	{
		alpha <- m * (m*(1-m)/v -1)
		beta <- (1-m) * (m*(1-m)/v -1) 
	}
	return(list(alpha=alpha,beta=beta))
}

##------------------------------addubg pseudocount
smoothCount <- function(freq.file,snp.file){
	d <- getSNP(freq.file, snp.file)
	d$Totalreadadj <- as.numeric(d$Totalread + 1)
	d$Altreadadj <- as.numeric( d$Altread + 1 )
	d$Altfreqadj <- as.numeric( d$Altreadadj / d$Totalreadadj) 
	return(d[,c("Chr","Pos","Totalreadadj","Altreadadj","Altfreqadj","Ale")])
}

getParaCase <- function(){
	##------------------------------empirical data
	return(list(alpha = 0.5,beta = 0.5))
}

#emission probability for model
getEmProb <- function(data,alpha,beta){
	## using the moment estimate
	# freq.adj <- smoothCount(n.file,snp.file)
	m <- apply(data,1,function(x){
				lchoose(x[1],x[2]) + lbeta(x[2]+alpha,x[1]-x[2]+beta)-lbeta(alpha,beta)
				#lbeta to tackle underflow
				# lgamma(x[1]+1) + lgamma(x[2]+alpha) + lgamma(x[1]-x[2]+beta)+lgamma(alpha+beta)
				# -lgamma(x[2]+1)-lgamma(x[1]-x[2]+1) -lgamma(x[1]+alpha+beta) - lgamma(alpha) - lgamma(beta)
				})
	return(m)
}

myVitebi <- function(tFreqFile,nFreqFile,snp.file){
	#test using 
	#observition space
	#oberservations: alt reads number, total read, and baf 
	tdata <- smoothCount(tFreqFile,snp.file)
	ndata <- smoothCount(nFreqFile,snp.file)
	tn <- merge(tdata,ndata,by=c("Chr","Pos"),suffixes=c("_t","_n"))
	tn <- tn[order(tn$Pos),]
	#state space 
	states <- c("N","D")
	kNumStates <-length(states)
	
	###------------------------------transition matrix
	tranMX <- matrix(c(1-1/50,1/50,1/20,1-1/20),
					nrow=2,ncol=2,byrow=TRUE,
					dimnames=list(states,states))

	##Emission Matrix------------------------------
	parmCase <- getParaCase()
	parmNull <- getParaNull(ndata)
	emiNomral <- getEmProb(tn[,c("Totalreadadj_n","Altreadadj_n")],parmNull$alpha,parmNull$beta)
	emiTumor <- getEmProb(tn[,c("Totalreadadj_t","Altreadadj_t")],parmCase$alpha,parmCase$beta)
	emiMX <- cbind("N"=emiNomral,"D"=emiTumor)
	row.names(emiMX) <- apply(tn,1,function(x){paste(x[1],x[2],sep="_")})
	
	#Initialization
	start <- c("N"=(1-1/50),"D"=1/50)
	num.windows <- dim(emiMX)[1]
	vtbM <- matrix(NA, nrow=num.windows, ncol=kNumStates)
	vtbPt <- matrix(NA, nrow=num.windows, ncol=kNumStates)
	vtbM[1,]<- log(start) + emiMX[1,]
	vtbPt[1,] <- c(0,0)
	##------------------------------Viterbi 
 	for(i in 2:num.windows){
 		tempI <- vtbM[(i-1),] + log(tranMX)
 		# if(i >= 134 && i <= 153){
 			# print(tempI)
 		# }
 		vtbM[i,]<- apply(tempI,1,max)
		vtbM[i,]<- vtbM[i,] + emiMX[i,]
		vtbPt[i,] <-  as.vector(apply(tempI,1,which.max))
	}	
	#output
	vtbStates <- matrix(NA, ncol=3,nrow=num.windows)
	vtbStates[,3] <- apply(vtbM,1,which.max)
	vtbStates[,2] <- tn$Pos
	vtbStates[,1] <- tn$Chr
	colnames(vtbStates) <- c("Chr","Pos","State")
  	return(vtbStates)
   }

plotVbtStates <- function(vtbStates,delReg) {
	reg <- read.delim(delReg)
	pdelR <- reg[which(reg$CHR == unique(vtbStates[,1])),] 
	# output <- paste(act.num,"_CNV_viterbi.pdf",sep="")
	require(ggplot2)
	# mytitle <- paste("CHR",unique(vtbStates[,1]),sep="")
	maxPos <- as.numeric(max(vtbStates$Pos))
  	plot(vtbStates$Pos,vtbStates$State,type="l"
  		, xlim=c(1,maxPos)
		,main=paste("CHR",unique(vtbStates[,1]),sep="")
		, xlab="",ylab="",xaxt="n",yaxt="n")
  	colors <- c("red","green","blue","pink","orange")
  	if(nrow(pdelR) > 1) {
		for (i in 1:nrow(pdelR)){
		abline(v=pdelR[i,3],col=colors[i])
		abline(v=pdelR[i,4],col=colors[i])
		}
  	}else if(nrow(pdelR) == 1){
	  	abline(v=pdelR[1,3],col=colors[1])
	  	abline(v=pdelR[1,4],col=colors[1])
	}
	
}

myVitebi22 <- function(act.num,acn.num){ 
	# act.num <- 13  # 10,13,16 
	# acn.num <- 19
	if(Sys.info()["sysname"] == "Darwin"){
		Sys.setenv(NETDIR="/Volumes/ys_lab_scratch/jh3283/net/"
				   ,SCRIPTS="/Users/jh3283/Dropbox/scripts/"
	      )
	}
	NETDIR <- Sys.getenv("NETDIR")
	SCRIPTS <- Sys.getenv("SCRIPTS")
	snp.file <- paste(NETDIR,"AC1SNPs.txt",sep="")
	CWD <- paste(NETDIR,"AC",act.num,"/res/",sep="")
	setwd(CWD)
	# source(paste(SCRIPTS,"projNET/testCNV/myBBViterbi.r",sep=""))
	freq_mode.tumor <- paste(NETDIR,"AC",act.num,"/BAF/","chrCHR.freq",sep="")
	freq_mode.normal <- paste(NETDIR,"AC",acn.num,"/BAF/","chrCHR.freq",sep="")

	##------------------------------test
	print(NETDIR)
	print(SCRIPTS)
	print(snp.file)
	print(freq_mode.tumor)
	print(freq_mode.normal)
	##------------------------------
	vtbStates <- data.frame(Chr=numeric(),Pos=numeric(),State=numeric())
	for (i in 1:22){
		tfile <- gsub("CHR",i,freq_mode.tumor)
		nfile <- gsub("CHR",i,freq_mode.normal)
		vtbStates <- rbind(vtbStates,myVitebi(tfile,nfile,snp.file))
	}
	print(table(vtbStates[,c(1,3)]))
	return(vtbStates)
}

plotVbtStates22 <- function(vtbStates,delReg.file,output){
	library(ggplot2)
	ppdata <- list()
	ncol <- 6
	nrow <- 4
	#test
	# output <- "test.pdf"
	delReg <- delReg.file
	pdf(output)
	par(mfrow=c(4,6),mai=c(0.1,0.3,0.3,0.1))
	for (i in 1:22) {
		print(plotVbtStates(vtbStates[vtbStates$Chr == i,],delReg))
  	}
  	mytitle <- gsub(".pdf","",output)
  	plot(1:10,type="n",xlab="",ylab="",xaxt="n",yaxt="n",axes=F)
  	text(5,5,unlist(strsplit(mytitle,"_"))[1],cex=1.5)
  	plot(1:10,type="n",xlab="",ylab="",xaxt="n",yaxt="n",axes=F)
  	text(5,5,unlist(strsplit(mytitle,"_"))[2],cex=1.5)
  	dev.off()

}




