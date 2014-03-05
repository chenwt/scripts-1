##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep. 17, 2013 
#Last Updated:
#Usage:
#Description: input a file contain gene names, one per line, out put uniprot accession, one per line
 # source("http://bioconductor.org/biocLite.R")
 # biocLite("biomaRt")

args <- commandArgs(TRUE)
 if (is.null(args)){
   print("Please provide parameters")
   exit
 }else{
   print(args)
 }
fn <- args[1]
# fn <- "/ifs/data/c2b2/ac_lab/jh3283/projMisc/yLan/Feb28/yLanGene.list"
gene <-unlist(read.table(fn,header=T))

print(paste("Total input gene name:",length(gene),sep=" "))

require(biomaRt)
mart <- useDataset(dataset="uniprot",mart=useMart("unimart"))
attributes <- c("accession","protein_name","gene_name","organism")
# attributes <- c("accession","gene_name")
uniprot <- getBM(filters="gene_name",mart=mart,values=gene,attributes=attributes)
uniprot <- uniprot[grep("Homo sapiens",uniprot$organism),1]
print(paste("Total output uniprot name: ",length(uniprot),sep=" "))

write.table(uniprot,paste(fn,"_unipro",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
