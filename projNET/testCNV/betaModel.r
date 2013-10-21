#!/bin/Rscirpts
#Author: Jing He
#Date: Mar.26,2013
#Last Update: Mar.26,2013
#Modelling the reads using beta distribution

##------------------------------null distribution
source("~/Dropbox/scripts/projNET/testCNV/CNVTest.r")

n.file <- gsub("CHR","10","/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chrCHR.freq")
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
	vtbStates <- vector(length=num.windows)
	vtbStates <- apply(vtbM,1,which.max)
	vtbStates$Pos <- tn$Pos
	vtbStates$Chr <- tn$Chr
  	return(vtbStates)
   }

delReg <- "/Volumes/ys_lab_scratch/jh3283/net/AC3somatic_corrected_NoModLong-3-type.seg"
plotVbtStates <- function(vtbStates,delReg) {
	delReg <- read.delim(delReg)
	pdelR <- delReg[which(delReg$CHR == vtbStates$Chr),] 
	output <- paste(act.num,"_CNV_viterbi.pdf",sep="")
	title <- paste("CHR",unique(vtbStates$Chr),sep="")
  	plot(vtbStates$Pos,vtbStates,type="l", 
  	main=title,xlab="Pos",ylab="States")
  	color <- c("red","green","blue","pink","orange")
  	for (i in 1:nrow(pdelR)){
	  	print(abline(v=pdelR[i,3],col=color[i]))
	  	print(abline(v=pdelR[i,4],col=color[i]))
	}

}

runViterb <- function(act.num,acn.num,snp.file,delReg.file){
	
  	tFreqFile <- gsub("CHR","13","/Volumes/ys_lab_scratch/jh3283/net/AC10/BAF/chrCHR.freq")
	nFreqFile <- gsub("CHR","13","/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chrCHR.freq")
	vtbStates <- myVitebi(tFreqFile,nFreqFile,snp.file)
}

