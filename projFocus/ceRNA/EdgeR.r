#!/usr/bin/Rscript
#Author: Jing He
#Date:24 Oct,2013 
#Last Updated:
#COMMENTS: need edgeR installed; 
#input: <string:path you wnat your results to be> 
# 		<string:name of your design file(4 cols, tab delimite:example)
#		<string:name of count matrix file>
# 		<string:name of your output files>
#output: <file:pdf of 2 plots> <file: txt of differetial expresssed genes>

# example files:
path <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp"
design_file <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/rnaSeq/input_edgeR_ColDesign.txt"
data_matrix_file <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/data/rnaSeq/brca_rnaseq_rawCount_436.mat"
output_file <- "EdgeR_run1"

argv <- commandArgs(TRUE)
if (length(argv) < 4) {
  cat("Usage: Rscript EdgeR.r path design_file data_matrix_file output_file number_of_replicate_for_condition_1 number_of_replicate_for_condition_2 ...\n")
  q(status = 1)
}

path <- argv[1]
design_file <- argv[2]
data_matrix_file <- argv[3]
output_file <- argv[4]

##----------------------------
#load data
dataMat <- read.delim(data_matrix_file,header=T)
rownames <- dataMat[,1]
dataMat <- dataMat[,-1]

datacolDesign <- read.table(design_file,sep= "\t")
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

designCol <- datacolDesign[,4]
names(designCol) <- datacolDesign[,1]
designCond <- na.omit(designCol[colnames(dataMat)])

dataMat <- dataMat[,names(designCond)]
condition <- factor( designCond )

##----------------------------
#fitting model
library(edgeR)
setwd(path)
y <- DGEList(counts=dataMat,group=condition)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
save(et,file=paste(output_file,".rda",sep=""))


###----------------------------
#output
pdf(paste(output_file,".pdf",sep=""))
plotBCV(y, cex=0.4)

de <- decideTestsDGE(et, p=0.05, adjust="BH")
detags <- rownames(topTags(et, n=500))
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")
dev.off()

outTag <- topTags(et,n=Inf,adjust.method="bonferroni",sort.by="p.value")
write.table(outTag,file=paste(output_file,"DEG.txt",sep=""),sep="\t",quote=F)
outTagMat <- dataMat[rownames(outTag),]
outTagMat_header <- as.matrix(t(c("Gene",colnames(outTagMat))))

write.table(outTagMat_header,file=paste(output_file,"DEGMat.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(outTagMat,file=paste(output_file,"DEGMat.txt",sep=""),sep="\t",quote=F,append = TRUE,col.names=F)
