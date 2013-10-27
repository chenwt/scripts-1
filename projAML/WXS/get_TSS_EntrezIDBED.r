#!/usr/bin/Rscript
#Author: Jing He
#Date:Sep. 17, 2013 
#Last Updated:
#Usage:
#input: <file: file first colum is the gene name> 
#output: <file: file column was replaced by entrezID, some rows might lost due to unmapping >
args <- commandArgs(TRUE)
 if (is.null(args)){
   print("Please provide parameters")
   exit
 }else{
   print(args)
 }
require(org.Hs.eg.db)
in_file <- args[1]
output_file <- paste("entrezID_",in_file,sep="")
# fn <- "input_test.txt"
data <- read.table(in_file,header=T,sep="\t")
gene <-unlist(data[,1])


geneSym <- org.Hs.egSYMBOL2EG
geneSym <- as.list(geneSym[mappedkeys(geneSym)])
entrez <- unlist(geneSym[gene])

data <- data[!is.na(entrez[gene]),]
data[,1] <- entrez[gene[!is.na(entrez[gene])]]
colnames(data) <- c("EntrezID",colnames(data)[-1])

write.table(data,output_file,row.names=F,quote=F,sep="\t")
