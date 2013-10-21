#!/bin/Rscirpts
## parameters include ac.num, snp.file "AC1SNPs_2.txt",
# cat AC1SNPs.txt | awk '{print $1"\t"$2+1"\t"$3}' > AC1SNPs_3.txt 
#example: Rscript ~/scripts/projNET/run_PlotFreq.r 13 AC1SNPs.txt 
#

args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters: ac#,and patient snp file name")
	exit
}else{
	ac.num <- args[1]
	snp.file <- args[2]
	print(c(ac.num, snp.file))
}

SCRIPTSDIR <- Sys.getenv("SCRIPTS")

###------------------------------ preparing filenames
NETDIR <- Sys.getenv("NETDIR")
print(NETDIR)
if(is.null(NETDIR))
{exit}

cwd <- paste(NETDIR,"AC",ac.num,"/res/",sep="")
print(cwd)
FreqFileMode <- paste(NETDIR,"AC",ac.num,"/BAF/chrCHR.freq",sep="")
print(c("Freq file MODE:", FreqFileMode))
output <- paste("AC",ac.num,"_22chr",
				as.character(unlist(strsplit(snp.file,"\\."))[[1]]),
				"_BAF.pdf",sep="")
print(c("output filename:", output))
snp.file <- paste(NETDIR,snp.file,sep="")
print(c("patient snp file:",snp.file))

###------------------------------ loading packages and source code, run functions
setwd(cwd)
source(paste(SCRIPTSDIR,"projNET/PlotFreq.r",sep=""))
require(ggplot2)

FreqBAFPlot22Chrs(FreqFileMode, output.file=output, snp.file, plot.title=paste("AC",ac.num,"BAF"))


