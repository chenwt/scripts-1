#!/bin/Rscirpts
#Author: Jing He
#Date: Mar.26,2013
#Last Update: Apr.15,2013
#Modelling the reads using beta distribution

##------------------------------null distribution
source("~/Dropbox/scripts/projNET/testCNV/CNVTest.r")

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
	d$Totalreadadj <- as.numeric(d$Totalread + 1)
	d$Altreadadj <- as.numeric( d$Altread + 1 )
	d$Altfreqadj <- as.numeric( d$Altreadadj / d$Totalreadadj) 
	return(d[,c("Chr","Pos","Totalreadadj","Altreadadj","Altfreqadj","Ale")])
}

##------------------------------get alpha & beta for null distribution
	##------------------------------
	# get alpha, beta value based on Chromosome level information
	# alpha, beta ~ AF ~ beta(alpha, beta)
	# need to integrate information from based level to zoom out expression nosiy.
getParaNull <- function(freq.adj){
	# get alphaBar, betaBar using getParaNullChr(Altfreqadj)
	Totalreadadj <- as.numeric(freq.adj[,1])
	Altfreqadj <- as.numeric(freq.adj[,2])
	aBar <- getParaNullChr(Altfreqadj)$alpha
	bBar <- getParaNullChr(Altfreqadj)$beta
	T <- 20
	alpha <- vector(length=length(Totalreadadj))
	beta <- vector(length=length(Totalreadadj))
	alpha <- unlist(lapply(Totalreadadj,function(x){
					min( sqrt(aBar * x / 2), T)}))
	beta <- unlist(lapply(Totalreadadj,function(x){
					min( sqrt(bBar * x / 2), T)}))
	return(list(alpha=alpha,beta=beta))
}

getParaNullChr <- function(af){
	m <- mean(af)
	v <- var(af)
	if(v < m * (1-m))
	{
		alpha <- m * (m * (1-m)/v -1)
		beta <- (1-m) * (m * (1-m)/v -1) 
	}
	return(list(alpha=alpha,beta=beta))
}

##--------------------------get alpha & beta for case distribution
getParaCase <- function(){
	return(list(alpha = 0.5,beta = 0.5))
}

#------------emission probability for model
getEmProbNull <- function(total,alt,baf,ap,bp){
	# Totalreadadj <- ReadsFile[,1]	
	# Altreadadj <- ReadsFile[,2]

	em <- apply(cbind(rep(a,length(ap)),rep(b,length(bp)),ap,bp),
			1,function(x){
						p <- getPai(x[1],x[2],x[3],x[4])
						print(p)
						log(myLikeliHood(p, x[1],x[2],x[3],x[4]))
						})

	return(em)
}


getEmProbDel <- function(totalt, altt, baft,totaln, altn, bafn,a,b,ap,bp){
	
	# Totalreadadj <- ReadsFile[,1]	
	# Altreadadj <- ReadsFile[,2]
	data <- cbind(rep(a,length(ap)),rep(b,length(bp)),ap,bp)
	em <- apply(data,1,function(x){
						p <- getPai(x[1],x[2],x[3],x[4])
						print(p)
						log(myLikeliHood(p, x[1],x[2],x[3],x[4]))
						})

	return(em)
}

getPai <- function(a,b,ap,bp){
	pai <- seq(0.3,1,by=0.1)
	p <- pai[which.max(unlist(lapply(pai,function(x)
			{myLikeliHood(x,a,b,ap,bp) })))]
	return(p)
}

myLikeliHood <- function(pai,a,b,ap,bp){
	# nm1 <- ceiling(pai * N)
	# mm1 <- ceiling(pai * m)
	# nm2 <-ceiling((1-pai) * N)
	# mm2 <-ceiling((1-pai) * m)
	lk <- pai * beta(a, b) + (1-pai) * beta(ap,bp)
		lchoose()
	return(lk)
}



myVitebi <- function(tFreqFile,nFreqFile,SNPFILE){
	#test using 
	#observition space
	#oberservations: alt reads number, total read, and baf 
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
	ap <- getParaNull(tn[,c("Totalreadadj_n","Altfreqadj_n")])$alpha
	bp <- getParaNull(tn[,c("Totalreadadj_n","Altfreqadj_n")])$beta
	a <- getParaCase()$alpha
	b <- getParaCase()$beta
	# parmCase <- getParaCase()
	# parmNull <- getParaNull(ndata)
	emiNomral <- getEmProb(a,b,ap,bp)
	
	ap2 <- getParaNull(tn[,c("Totalreadadj_t","Altfreqadj_t")])$alpha
	bp2 <- getParaNull(tn[,c("Totalreadadj_t","Altfreqadj_t")])$beta
	emiTumor <- getEmProb(a,b,ap2,bp2)
	
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




