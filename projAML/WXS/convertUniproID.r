##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep. 17, 2013 
#Last Updated:
#Usage:
#Description: get the genename of uniprob id of proteins in preppi database

require(biomaRt)
mart <- useDataset(dataset="uniprot",mart=useMart("unimart"))
data <-read.table("~/SCRATCH/database/preppi/preppi_int_3col.txt",header=T)

NodeUnipro <- c(unique(unlist(data$int_a)),unique(unlist(data$int_b)))
unipro2genename <- getBM(filters="accession",mart=unipro,values=NodeUnipro,attributes=c("accession","gene_name"))
write.table(unipro2genename,"~/SCRATCH/database/geneIDconverter/unipro2genenamePrePPI.txt",row.names=F,quote=F,sep="\t")
save(unipro2genename,file="~/SCRATCH/database/geneIDconverter/unipro2genenamePrePPI.rda")