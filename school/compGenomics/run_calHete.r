#!/bin/Rscirpts
#Author: Jing He
#Date: Mar. 2013
#Last Update: Mar.24, 2013
#parameters: ac.num snp.file
#example: Rscript ~/scripts/projNET/run_calHete.r 13 AC1SNPs.txt  

## cat AC1SNPs.txt | awk '{print $1"\t"$2+1"\t"$3}' > AC1SNPs_2.txt 

##------------------------------header
args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters: ac#,and patient snp file name")
	exit
}else{
	ac.num <- args[1]
	snp.file <- args[2]
	print(c(ac.num, snp.file))
}

Sys.setenv(	NETDIR="/Volumes/ys_lab_scratch/jh3283/net/",
			SCRIPTS="/Users/jh3283/Dropbox/scripts/"
			)



SCRIPTSDIR <- Sys.getenv("SCRIPTS")
source(paste(SCRIPTSDIR,"projNET/BAFTest.r",sep=""))
###------------------------------testing on local
# ac.num <- 13
# NETDIR <- "/Volumes/ys_lab_scratch/jh3283/net/"
# SCRIPTSDIR <- "/Users/jh3283/Dropbox/scripts/"
# snp.file <- "AC1SNPs.txt"
###------------------------------ preparing filenames
NETDIR <- Sys.getenv("NETDIR")
print(NETDIR)
if(is.null(NETDIR))
{exit}

cwd <- paste(NETDIR,"AC",ac.num,"/res/",sep="")
print(cwd)
FreqFileMode <- paste(NETDIR,"AC",ac.num,"/BAF/chrCHR.freq",sep="")
# print(c("Freq file MODE:", FreqFileMode))
outname <- paste("CNV_AC",ac.num,"_22chr_",
				as.character(unlist(strsplit(unlist(strsplit(snp.file,"/"))[length(unlist(strsplit(snp.file,"/")))],"\\."))[[1]]),
				".pdf",sep="")
output <- paste(cwd,outname,sep="")
# print(c("output filename:", output))
snp.file <- paste(NETDIR,snp.file,sep="")
# print(c("patient snp file:",snp.file))
# plotCNV22("/Volumes/ys_lab_scratch/jh3283/net/AC13/BAF/chrCHR.freq","/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt",paste("CNV_",ac.num,".pdf",sep=""))
# plotCNV22("/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chrCHR.freq","/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt",paste("CNV_",ac.num,".pdf",sep=""))

plotCNV22(FreqFileMode,snp.file,output)


# snp.file <-  "/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt"
# t.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC13/BAF/chr11.freq"
# n.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chr11.freq"
# bin.num <- 50
# install.packages("reshape2")
# require(reshape2)
# dcast(SNPdata, cut(AltFreq, breaks = c(0, 25, 100)) ~ status)

# x <- as.numeric(SNPdata$Altfreq )
# br <- seq(0, max(SNPdata$Altfreq), by=0.05)
# xt <- summary(cut(x[x!=0], br, labels=br[-1], include.lowest=T, ordered_result=T))

# vapply(SNPdata$Altfreq,mySNPHete,1)
# library(help="diptest")
# dip.test(SNPdata$Altfreq[SNPdata$Altfreq > 0])
# plot(bkde(v[61:70]),type="l",cex=3)



