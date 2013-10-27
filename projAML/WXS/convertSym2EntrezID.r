#!/usr/bin/Rscript
#Author: Jing He
#Date:Sep. 17, 2013 
#Last Updated:
#Usage:
#input: <file: gene name by line> 
#output: <Entrez ID, by line >
args <- commandArgs(TRUE)
 if (is.null(args)){
   print("Please provide parameters")
   exit
 }else{
   print(args)
 }
require(org.Hs.eg.db)
in_file <- args[1]
# fn <- "input_test.txt"
output_file <- paste("mappped.",in_file,sep="")
gene <-unlist(read.table(in_file,header=T))

geneSym <- org.Hs.egSYMBOL2EG
geneSym <- as.list(geneSym[mappedkeys(geneSym)])
out <- cbind(gene, unlist(geneSym[gene]))
colnames(out) <- c("GeneName","EntrezID")
# print(paste("input number :",length(gene),"output number:", dim(out)[1],sep="    "))

write.table(out,output_file,row.names=F,quote=F,sep="\t")
