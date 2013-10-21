#!/bin/Rscirpts
#Author: Jing He
#Date: Mar.26,2013
#Last Update: Apr.15,2013
#Modelling the reads using beta binomial distribution

##------------------------------null distribution
# source("~/Dropbox/scripts/projNET/testCNV/CNVTest.r")

getSNP <- function(freq.file,snp.file){
  freq <- read.table(freq.file, skip=5, header=TRUE)
  freq$Altfreq <- freq$Altreads / freq$Totalreads
  freq <- freq[freq$Altfreq > -1 & freq$Altfreq < 2, ]
  snp <- read.table(snp.file,sep="\t",header=F)
  colnames(snp) <- c("Chr","Pos","Ale")
  #------------------------------ load 22 freq data
  snpFreq <- merge(x=freq,y=snp,by=c("Chr","Pos"))
  snpFreq$Altfreq <- snpFreq$Altreads / snpFreq$Totalreads
  return(snpFreq)
}

##------------------------------reading data and adding pseudocount
smoothCount <- function(freq.file,SNPFILE){
	d <- getSNP(freq.file, SNPFILE)
	d$Totalreadadj <- as.numeric(d$Totalread + 2)
	d$Altreadadj <- as.numeric( d$Altread + 1 )
	d$Altfreqadj <- as.numeric( d$Altreadadj / d$Totalreadadj) 
	return(d[,c("Chr","Pos","Totalreadadj","Altreadadj","Altfreqadj")])
}

getParaNull <- function(freq.adj){
	require(MASS)

	Totalreadadj <- as.numeric(freq.adj[,1])
	Altfreqadj <- as.numeric(freq.adj[,2])
	aBar <- getParaNullChr(Altfreqadj)$alpha
	bBar <- getParaNullChr(Altfreqadj)$beta
	T <- 20
	alpha <- vector(length=length(Totalreadadj))
	beta <- vector(length=length(Totalreadadj))
	alpha <- sapply(Totalreadadj,function(x){
					min( sqrt(aBar * x / 2), T)})
	beta <- unlist(lapply(Totalreadadj,function(x){
					min( sqrt(bBar * x / 2), T)}))
	return(list(alpha=alpha,beta=beta))
}

getParaNullChr <- function(af){
	m <- mean(af)
	v <- var(af)
	if(v < m * (1-m))
	{
		alpha <- m * (m*(1-m)/v -1)
		beta <- (1-m) * (m*(1-m)/v -1) 
	}
	return(list(alpha=alpha,beta=beta))
}

getParaDel <- function(paramNull){
	a0 <- sapply(paramNull$alpha,as.numeric)
	b0 <- sapply(paramNull$beta,as.numeric)

	a1 <- sapply(1:length(a0),
			FUN=function(x){
				mean(sapply(seq(0.001,.999,by=0.001),
					FUN=function(y){(1-y)/(2-y) * (a0[x] + b0[x]) * y}), 
					0.5) })
	b1 <- sapply(1:length(a0),
		function(x){a0[x] + b0[x] - a1[x]})

	return(list(alpha1 = a1,beta1 = b1,alpha2=b1, beta2=a1))
}

getEmProbDel <- function(baf,paramDel){
	baft <- as.numeric(baf)
	a1 <- as.numeric(paramDel$alpha1)
	a2 <- as.numeric(paramDel$alpha2)
	b1 <- as.numeric(paramDel$beta1)
	b2 <- as.numeric(paramDel$beta2)
	em <-sapply(1:length(a1),function(i){
					0.5 * dbeta(baft[i],a1[i],b1[i]) +
					0.5 * dbeta(baft[i],a2[i],b2[i])
					})
	
	return(log(em))		
}

getEmProbNull <- function(baf,ap,bp) apply(cbind(baf,ap,bp),1,function(x){log(dbeta(x[1],x[2],x[3]))})


myVitebi <- function(tFreqFile,nFreqFile,SNPFILE){
	#observation
	tdata <- smoothCount(tFreqFile,SNPFILE)
	ndata <- smoothCount(nFreqFile,SNPFILE)
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
	
    paraNormal <- getParaNull(tn[,c(6,8)])
    paraDel <- getParaDel(paraNormal)
    emiTumor <- getEmProbDel(tn$Altfreqadj_t,paraDel)
    emiNormal <- getEmProbNull(tn$Altfreqadj_n,paraNormal$alpha,paraNormal$beta)
	emiMX <- cbind("N"=emiNormal,"D"=emiTumor)
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
  	return(as.data.frame(vtbStates))
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
	SNPFILE <- paste(NETDIR,"AC1SNPs.txt",sep="")
	CWD <- paste(NETDIR,"AC",act.num,"/res/",sep="")
	setwd(CWD)
	# source(paste(SCRIPTS,"projNET/testCNV/myBBViterbi.r",sep=""))
	freq_mode.tumor <- paste(NETDIR,"AC",act.num,"/BAF/","chrCHR.freq",sep="")
	freq_mode.normal <- paste(NETDIR,"AC",acn.num,"/BAF/","chrCHR.freq",sep="")

	##------------------------------test
	print(NETDIR)
	print(SCRIPTS)
	print(SNPFILE)
	print(freq_mode.tumor)
	print(freq_mode.normal)
	##------------------------------
	vtbStates <- data.frame(Chr=numeric(),Pos=numeric(),State=numeric())
	for (i in 1:22){
		tfile <- gsub("CHR",i,freq_mode.tumor)
		nfile <- gsub("CHR",i,freq_mode.normal)
		vtbStates <- rbind(vtbStates,myVitebi(tfile,nfile,SNPFILE))
	}
	# print(table(vtbStates[,c(1,3)]))
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



