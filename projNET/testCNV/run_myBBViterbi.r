#!/bin/Rscirpts
#Author: Jing He
#Date: Mar 28, 2013
#Last Update: Mar 28, 2013
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
##------------------------------test 
act.num <- 13
acn.num <- 19

if(Sys.info()["sysname"] == "Darwin"){
		Sys.setenv(NETDIR="/Volumes/ys_lab_scratch/jh3283/net/"
				   ,SCRIPTS="/Users/jh3283/Dropbox/scripts/"
	      )
	}
NETDIR <- Sys.getenv("NETDIR")
SCRIPTS <- Sys.getenv("SCRIPTS")
source(paste(SCRIPTS,"projNET/testCNV/myBBViterbi.r",sep=""))

runViterb <- function(act.num, acn.num){
	# FREQFILEMODE <- 
	states <- myVitebi22(act.num,acn.num)
	delReg <- "/Volumes/ys_lab_scratch/jh3283/net/AC3somatic_corrected_NoModLong-3-type.seg"
	# plotVbtStates(vtbStates,delReg)
	out <- paste("AC",act.num,"_vtbCNV.pdf",sep="")
	plotVbtStates22(states, delReg, out)
}




