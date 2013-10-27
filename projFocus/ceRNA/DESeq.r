#!/usr/bin/Rscript
#J.HE
#input:
#output:
#TODO:

 # args <- commandArgs(TRUE)
 # if (is.null(args)){
 #   print("Please provide parameters")
 #   exit
 # }else{
 #   print(args)
 # }

# library(EBSeq)
library(DESeq)
rootd <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp"
setwd(rootd)
# load raw expression matrix and prepare a matrix
# infile <- args[1]
infile <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/rnaSeq2/brca_rnaSeq2_rsem_raw_1.mat"
designfile <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/rnaSeq2/input_EBseqR_ColDesign.txt"
datacolDesign <- read.table(designfile,	sep= "\t")

# 
table(datacolDesign[,4])
replace = function(df,old,new){
	dfnew <- apply(df,c(1,2),function(x){
	pattern <- paste("^",old,"$",sep="")
	gsub(pattern,new,gsub(" ","",x),perl=T)
	})
	return(dfnew)
}
# take tumor and normal sample
datacolDesign <- datacolDesign[which(datacolDesign[,4]=="1" | datacolDesign[,4]=="11"),]
datacolDesign <- replace(datacolDesign,1,"tumor")
datacolDesign <- replace(datacolDesign,11,"normal")
datacolDesign[,1] <- sapply(datacolDesign[,1],function(x){gsub("-",".",x)})


dataMat <- read.table(infile,sep="\t",header=T)
row.names(dataMat) <- dataMat[,1]
dataMat <- dataMat[,-1]

designCol <- datacolDesign[,4]
names(designCol) <- datacolDesign[,1]
designCond <- na.omit(designCol[colnames(dataMat)])
dataMat <- dataMat[,names(designCond)]

# prepare the Design
libDesign <- rep("paired-end",ncol(dataMat))
dataMatDesign = data.frame(row.names = colnames(dataMat),
 		condition = designCond,
		libType = libDesign )
condition <- factor( designCond )

# normalization DESeq
dataMat <- apply(dataMat,c(1,2),round)
cds = newCountDataSet( dataMat, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
topleft( counts( cds, normalized=TRUE ) )
cds = estimateDispersions( cds )
str( fitInfo(cds) )
save(cds, file=paste(rootd,"/DESeq_run1_cds.rda",sep=""))

head( fData(cds) )
res = nbinomTest( cds, "tumor", "normal" )
save(res, file=paste(rootd,"/DESeq_run1_res.rda",sep=""))
head(res)



pdf(paste(rootd,"/DESeq_run1.pdf",sep=""))
 plotDispEsts( cds )
 plotMA(res)
 hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

# resSig = res[ res$padj < 0.1, ]
# head( resSig[ order(resSig$pval), ] )
# head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )
write.csv( res, file=paste(rootd,"/DESeq_run1_DEG_result.csv",sep=""))

# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")




