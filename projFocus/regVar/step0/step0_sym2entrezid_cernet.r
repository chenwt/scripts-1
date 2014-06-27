cernafile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_ceRNA_network.txt"

cernet <- read.table(cernafile, header=T,sep="\t", stringsAsFactors=F)


# biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

## Get the entrez gene identifiers that are mapped to a gene symbol (and viceversa
mapped_genes <- mappedkeys(org.Hs.egSYMBOL2EG)
list_symbol2eg <- as.character(org.Hs.egSYMBOL2EG[mapped_genes])
mapped_genes <- mappedkeys(org.Hs.egSYMBOL)
list_eg2symbol <- as.character(org.Hs.egSYMBOL[mapped_genes])


cernetNew <- cernet

cernetNew$RNA1 <- list_symbol2eg[cernet[,1]]
cernetNew$RNA2 <- list_symbol2eg[cernet[,2]]
head(cernetNew)
write.table(cernetNew,file=paste(cernafile,"_entrezid.txt",sep=""), quote=F, row.names=F,col.names=T,sep="\t")



###-----test
list_symbol2eg['AAK1']
