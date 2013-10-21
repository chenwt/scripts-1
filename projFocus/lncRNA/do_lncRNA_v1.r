


#TODO: run this code sucessfully!
# rootd <- "/ifs/scratch/c2b2/ac_lab/jh3283/"
rootd <- "/ifs/scratch/c2b2/ac_lab/jh3283/"

library(biomaRt)
# listMarts()

mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
id_lncrna <- unique(read.table(paste(rootd,"ref/gencode/gencode.v18.lncRNA_transcripts_geneid_genename.txt",sep=""),
				heade=T))
ensembl_genes<- unlist(id_lncrna$gene_id) ## etc
gene_name <- unlist(id_lncrna$gene_name)
# entrezgene <- c("100130426","343066")
genename2entrez <- getBM(
        attributes = c("hgnc_symbol", "entrezgene"), 
        filters = "hgnc_symbol",
	    values= gene_name,
	    mart= mart)

# genename2ensembl <- getBM(
#         attributes = c("hgnc_symbol", "ensembl_gene_id"), 
#         filters = "hgnc_symbol",
# 	    values= gene_name,
# 	    mart= mart)

genename2entrez <- na.omit(genename2entrez)



##----------------------------
#loading level3 data
datad <- "projFocus/lncRNA/tcga_prad_RNAseqV2_Level_3/"
setwd(paste(rootd,datad,sep=""))
fns <- dir(pattern="genes.results")
readtable <- function(f){
	tmp <- read.table(f,header=T,sep="\t")
	tmp$gene_id <- sapply(tmp$gene_id,
		function(x){unlist(strsplit(x,split="\\|"))[2]})
	row.names(tmp) <- tmp$gene_id
	return(tmp)
}

for (i in 1:length(fns)){
	tmp <- readtable(fns[i])
	tmp_lncrna <- tmp[genename2entrez$entrezgene,2]
	if(i == 1) {lncrna_mat <- tmp_lncrna
	}else{
		lncrna_mat <- cbind(lncrna_mat,tmp_lncrna)
	}
}
# study design information
barcode <- read.table(
	pipe("cut -f2,22 ~/SCRATCH/projFocus/lncRNA/tcga_prad_RNAseqV2_Level_3_clinical/METADATA/UNC__IlluminaHiSeq_RNASeqV2/unc.edu_PRAD.IlluminaHiSeq_RNASeqV2.1.7.0.sdrf.txt"),sep="\t")
barcode <- barcode$V1[which(fns %in% barcode$V2 ==T,arr.ind=T)]

colnames(lncrna_mat) <- sapply(fns,function(x) barcode$V1[grep(x,barcode$V2)])

save(lncrna_mat,file="lncRNA_prad_rawcount.rda")


##----------------------------
# normalize data

# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq")

design <- rep("NA",length(fns))
design[grep("(01A)|(01B)",colnames(lncrna_mat),perl=T)] <- 
	rep("tumor",length(grep("(01A)|(01B)",colnames(lncrna_mat),perl=T)))
design[grep("(11A)|(11B)",colnames(lncrna_mat),perl=T)] <- 
	rep("normal",length(grep("(11A)|(11B)",colnames(lncrna_mat),perl=T)))
design <- factor(design)
cds = newCountDataSet(lncrna_mat, design )


