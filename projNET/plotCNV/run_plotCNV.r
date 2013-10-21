#!/bin/Rscirpts
#Author: Jing He
#Date: Mar.24,2013
#Last Update: Mar.24,2013
#parameters: ac.num snp.file
#example: Rscript ~/scripts/projNET/run_plotCNV.r 13 AC1SNPs.txt  

## cat AC1SNPs.txt | awk '{print $1"\t"$2+1"\t"$3}' > AC1SNPs_2.txt 

##------------------------------header
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters: ac#,and somatics snp file name")
	exit
}else{
	ac.num <- args[1]
	snp.file <- args[2]
	print(c(ac.num, snp.file))
}

if(Sys.info()["sysname"] == "Darwin") {
	Sys.setenv(	NETDIR="/Volumes/ys_lab_scratch/jh3283/net/",
				SCRIPTS="/Users/jh3283/Dropbox/scripts/"
				)
}


NETDIR <- Sys.getenv("NETDIR")
# print(NETDIR)
# if(is.null(NETDIR))
# {exit}
SCRIPTSDIR <- Sys.getenv("SCRIPTS")
source(paste(SCRIPTSDIR,"projNET/plotCNV.r",sep=""))
###------------------------------testing on local
# ac.num <- 13
# NETDIR <- "/Volumes/ys_lab_scratch/jh3283/net/"
# SCRIPTSDIR <- "/Users/jh3283/Dropbox/scripts/"
# snp.file <- "AC1SNPs.txt"
###------------------------------ preparing filenames

cwd <- paste(NETDIR,"AC",ac.num,"/res/",sep="")
print(paste("Working directory: ",cwd,sep=""))

FreqFileMode <- paste(NETDIR,"AC",ac.num,"/BAF/chrCHR.freq",sep="")
# print(c("Freq file MODE:", FreqFileMode))
outname <- paste("CNV_AC",ac.num,"_22chr_",
				as.character(unlist(strsplit(unlist(strsplit(snp.file,"/"))[length(unlist(strsplit(snp.file,"/")))],"\\."))[[1]]),
				".pdf",sep="")
output <- paste(cwd,outname,sep="")
# print(c("output filename:", output))
snp.file <- paste(NETDIR,snp.file,sep="")
# print(c("patient snp file:",snp.file))
plotCNV22(FreqFileMode,snp.file,output)


