#!/bin/Rscirpts
#Author: Jing He
#Date: Mar 28, 2013
#Last Update: Apr. 9, 2013
#parameters: 
#example: Rscript ~/scripts/projNET/testCNV/run_myBBViterbi.r 13 19


###------------------------------header------------------------------
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters")
	exit
}else{
	print(args)
}

###------------------------------start coding##------------------------------
act.num <- args[1]
acn.num <- args[2]


if(Sys.info()["sysname"] == "Darwin"){
		Sys.setenv(NETDIR="/Volumes/ys_lab_scratch/jh3283/net/"
				   ,SCRIPTS="/Users/jh3283/Dropbox/scripts/"
	      )
	}
NETDIR <- Sys.getenv("NETDIR")
SCRIPTS <- Sys.getenv("SCRIPTS")
SNPFILE <- paste(NETDIR,"AC1SNPs.txt",sep="")
DELREG <- paste(NETDIR,"AC3somatic_corrected_NoModLong-3-type.seg",sep="")

source(paste(SCRIPTS,"projNET/testCNV/myBBViterbi_v1.r",sep=""))

##------------------------------test 
ac.num <- 13
acn.num <- 19
CHR <- "CHR"
CHR <- 10
tFREQFILEMODE <- paste(NETDIR,"AC",ac.num,"/BAF/chr",CHR,".freq",sep="")
nFREQFILEMODE <- paste(NETDIR,"AC",acn.num,"/BAF/chr",CHR,".freq",sep="")

#chromosome level
#vtbStates <- myVitebi(tFREQFILEMODE,nFREQFILEMODE,SNPFILE)
#plotVbtStates(vtbStates,DELREG)

vtbsts<- myVitebi22(13,19)
plotVbtStates22(vtbsts,DELREG,"AC13_V21.pdf")
vtbsts<- myVitebi22(10,19)
plotVbtStates22(vtbsts,DELREG,"AC10_V21.pdf")
vtbsts<- myVitebi22(16,19)
plotVbtStates22(vtbsts,DELREG,"AC16_V21.pdf")
vtbsts<- myVitebi22(19,19)
plotVbtStates22(vtbsts,DELREG,"AC19_V21.pdf")

##------------------------------testing
#tFreqFile <- tFREQFILEMODE
#nFreqFile <- nFREQFILEMODE
#tdata <- smoothCount(tFreqFile,SNPFILE)
#ndata <- smoothCount(nFreqFile,SNPFILE)
#tn <- merge(tdata,ndata,by=c("Chr","Pos"),suffixes=c("_t","_n"))
#tn <- tn[order(tn$Pos),]
#
#paraNormal <- getParaNull(tn[,c(6,8)])
#paraDel <- getParaDel(paraNormal)
#emiDel <- getEmProbDel(tn$Altfreqadj_t,paraDel)
#emiNull <- getEmProbNull(tn$Altfreqadj_n,paraNormal$alpha,paraNormal$beta)
#
#
#nFreqFileltfreqadj_t,paraDel)
#dataSample <- getSnpFreq22Data(FREQFILEMODE,SNPFILE)
#head(dataSample)
#dim(dataSample)
#
#CHR <- 13
#FREQFILEMODE <- paste(NETDIR,"AC",ac.num,"/BAF/chr",CHR,".freq",sep="")
#dataChr <- smoothCount(FREQFILEMODE, SNPFILE)
#param <- getParaNull(dataChr[,3:4])$alpha
#summary(param)
#
#
#CHR <- 13
#FREQFILEMODE <- paste(NETDIR,"AC",acn.num,"/BAF/chr",CHR,".freq",sep="")
#dataChr1 <- smoothCount(FREQFILEMODE, SNPFILE)
#param1 <- getParaNull(dataChr1[,3:4])$alpha
#summary(param1)
#plot(1:nrow(dataChr),param,type="l")
#lines(param1,col="blue")
#getEmProbDel(cbind(dataChr$Altfreqadj,dataChr1$Altfreqadj),param,param,param1,param1)
#
#
#m <- mean(tdata$Altfreqadj)
#v <- var(tdata$Altfreqadj)
#alpha <- m * (m*(1-m)/v -1)
#beta <- alpha * (1/m-1)
#plot(density(rbeta(100,alpha,beta)))
#
###------------------------------plot reads count <==> freq 
#
#plot(dataSample$Altreads,dataSample$Altfreq,pch=20)
#title(main="AC13")
#
#runViterb <- function(act.num, acn.num){
#	# FREQFILEMODE <- 
#	states <- myVitebi22(act.num,acn.num)
#	delReg <- "/Volumes/ys_lab_scratch/jh3283/net/AC3somatic_corrected_NoModLong-3-type.seg"
#	# plotVbtStates(vtbStates,delReg)
#	out <- paste("AC",act.num,"_vtbCNV.pdf",sep="")
#	plotVbtStates22(states, delReg, out)
#}
#
#
#
#
#testEmissionDel <- function(tFreqFile,nFreqFile,SNPFILE){
#	tdata <- smoothCount(tFreqFile,SNPFILE)
#	ndata <- smoothCount(nFreqFile,SNPFILE)
#	tn <- merge(tdata,ndata,by=c("Chr","Pos"),suffixes=c("_t","_n"))
#	tn <- tn[order(tn$Pos),]
#	##Emission Matrix------------------------------
#	ap <- getParaNull(tn[,c("Totalreadadj_n","Altfreqadj_n")])$alpha
#	bp <- getParaNull(tn[,c("Totalreadadj_n","Altfreqadj_n")])$beta
#	a <- getParaCase()$alpha
#	b <- getParaCase()$beta
#	emiNomral <- getEmProbNull(tn$Altfreqadj_n,ap,bp)
#	emiTumor <- getEmProbDel(tn[,c("Altfreqadj_t","Altfreqadj_n")],a,b,ap,bp)
#	emiMX <- cbind("N"=emiNomral,"D"=emiTumor)
#	row.names(emiMX) <- apply(tn,1,function(x){paste(x[1],x[2],sep="_")})
#	return(emiMX)
#}
