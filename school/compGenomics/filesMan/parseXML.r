#!/bin/Rscirpts
#Author: Jing
#Date:  03.13.2013
#parameters: 
#example: Rscript ~parseXML.r file.xml

###------------------------------header------------------------------
arg <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters: file.xml")
	exit
}else{
	print(arg)
}

###------------------------------start coding##------------------------------
require("XML")
# print(arg)
disea_abb <- as.character(strsplit(arg,"\\.")[[1]])[1] 
# wdir <- "" 
wdir <- "/Volumes/ys_lab_scratch/jh3283/school/compGenomic/cgquery/" 
setwd(wdir)
 disea_abb <- "GBM_WXS"
# print(disea_abb)
# print(paste(wdir,disea_abb,".xml",sep=""))
fname <- paste(wdir,disea_abb,".xml",sep="")
dataDF <- xmlToDataFrame(fname) 
doc <- xmlParse(fname,useInternalNodes=TRUE)
doci <- xmlInternalTreeParse(fname)
filesizeLS <- xpathApply(doc,"//filesize")
filenameLS <- xpathApply(doc,"//filename")
sample_type <- xpathApply(doc,"//sample_type")
analysis_id <- xpathApply(doc,"//analysis_id")
refassem_short_name <- xpathApply(doc,"//refassem_short_name")
resDF <- t(sapply(1:length(filenameLS),function(x){cbind(xmlValue(analysis_id[[x]]), xmlValue(filenameLS[[x]]),
                                          xmlValue(filesizeLS[[x]]), xmlValue(sample_type[[x]]),
                                          xmlValue(refassem_short_name[[x]]) )}))
colnames(resDF) <- c("analysis_id","filename","filesize","sample_type","refassem_short_name")


write.table(resDF[,1],paste(wdir,disea_abb,"_analysis_id.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(resDF,paste(wdir,disea_abb,".txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
