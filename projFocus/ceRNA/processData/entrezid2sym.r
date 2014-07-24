# biocLite("org.Hs.eg.db")
rm(list=ls())
library(org.Hs.eg.db)
file = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_geneid-regulon.rda.4col.txt"
data = read.table(file,sep="\t",stringsAsFactors=F, header =T,colClasses="character")

## Get the entrez gene identifiers that are mapped to a gene symbol (and viceversa
mapped_genes <- mappedkeys(org.Hs.egSYMBOL2EG)
list_symbol2eg <- as.character(org.Hs.egSYMBOL2EG[mapped_genes])
mapped_genes <- mappedkeys(org.Hs.egSYMBOL)
list_eg2symbol <- as.character(org.Hs.egSYMBOL[mapped_genes])


dataold <-data
data$tf <- list_eg2symbol[data[,1]]
data$tfmod <- list_eg2symbol[data[,2]]

##check
setdiff(which(data$tf=='GATA3'), which(dataold$tf=='2625'))

write.table(na.omit(data),file=paste(file,".naomit.symbol",sep=""), quote=F, row.names=F,col.names=T,sep="\t")
