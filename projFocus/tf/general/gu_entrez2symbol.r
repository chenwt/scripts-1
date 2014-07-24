rm(list=ls())
file = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_geneid-regulon.rda.4col.txt"

# file = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/data/networks/brca_tcga_rnaseq851_geneid-regulon.rda"
# load(file)

dataDF = read.table(file,sep="\t", stringsAsFactors=F,header=T)
library(org.Hs.eg.db)
mapped_genes <- mappedkeys(org.Hs.egSYMBOL2EG)
list_symbol2eg <- as.character(org.Hs.egSYMBOL2EG[mapped_genes])
mapped_genes <- mappedkeys(org.Hs.egSYMBOL)
list_eg2symbol <- as.character(org.Hs.egSYMBOL[mapped_genes])

tfcol1  = dataDF$tf
tfcol2  = dataDF$tfmod

tfcol1_symb = list_eg2symbol[tfcol1]
tfcol2_symb = list_eg2symbol[tfcol2]
dataDF$tf <- tfcol1_symb
dataDF$tfmod <- tfcol2_symb


outfile = paste(file, ".symbol",sep="")
write.table(dataDF, outfile, sep="\t", col.names=T,row.names=F,quote=F)

outfile = paste(file, ".naomit.symbol",sep="")
write.table(na.omit(dataDF), outfile, sep="\t", col.names=T,row.names=F,quote=F)
