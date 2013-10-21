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
wdir <- "" 
# wdir <- "/Volumes/ys_lab_scratch/jh3283/school/compGenomic/" 
# disea_abb <- "GBM_exome"
# print(disea_abb)
# print(paste(wdir,disea_abb,".xml",sep=""))
laml.df <- xmlToDataFrame(paste(wdir,disea_abb,".xml",sep="")) 
laml.df <- laml.df[3:nrow(laml.df),] 
laml.df <- laml.df[!is.na(laml.df$analysis_id),]  
write.table(laml.df$analysis_id,paste(wdir,disea_abb,"_analysis_id.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(laml.df,paste(wdir,disea_abb,".txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
