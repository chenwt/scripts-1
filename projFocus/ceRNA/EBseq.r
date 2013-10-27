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

library(EBSeq)
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

# normalization EBSeq
Sizes = MedianNorm(dataMat)
dataMat = data.matrix(dataMat)
EBOut = EBTest(Data=dataMat,
 		Conditions=condition,
 		sizeFactors=Sizes, 
 		maxround=5,
 		Pool=TRUE)
PP=GetPPMat(EBOut)
DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]

DEfoundMat=dataMat[DEfound,]
save(EBOut, file=paste(rootd,"/EBSeq_run1_EBout.rda",sep=""))
save(DEfoundMat, file=paste(rootd,"/EBSeq_DEfoundMat.rda",sep=""))


GeneFC=PostFC(EBOut)
pdf(paste(rootd,"/EBSeq_run1.pdf",sep=""))
PlotPostVsRawFC(EBOut,GeneFC)
par(mfrow=c(2,2))
QQP(EBOut)
par(mfrow=c(2,2))
DenNHist(EBOut)
dev.off()

# resSig = res[ res$padj < 0.1, ]
# head( resSig[ order(resSig$pval), ] )
# head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] )
write.csv(DEfoundMat, file=paste(rootd,"/DESeq_run1_DEG_result.csv",sep=""))

# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")




