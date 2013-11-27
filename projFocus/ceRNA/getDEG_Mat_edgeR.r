#!/usr/bin/Rscript
#Author: Jing He
#Date:24 Oct,2013 
#Last Updated:
#COMMENTS: need edgeR installed; 
#input: <string:path you wnat your results to be> 
# 		<string:name of your design file(4 cols, tab delimite:example)
#		<string:name of count matrix file>
# 		<string:name of your output files>
#output: <file:pdf of 3 plots> <file: txt of differetial expresssed genes> <file: mat file of normalized row counts>

# # -------------when running local-------
# path <- "/Volumes/ac_lab/jh3283/SCRATCH/projFocus/ceRNA/result/exp/"
# setwd(path)
# design_file <- paste(path,"input_edgeR_ColDesign.txt",sep="")
# output_file <- "edgeR_brca_l3"
# data_matrix_file <- paste(path,"brca_rnaseq_rawCount_l3.mat",sep="")
# #---------------end_with_running_local-------

#argv <- commandArgs(TRUE)
#if (length(argv) < 4) {
#  cat("Usage: Rscript ~/scripts/projFocus/ceRNA/getDEG_Mat_edgeR.r  <working_path> <design_file> <data_matrix_file> <output_file> number_of_replicate_for_condition_1 number_of_replicate_for_condition_2 ...\n")
#  q(status = 1)
#}
#
#path <- argv[1]
#design_file <- argv[2]
#data_matrix_file <- argv[3]
#output_file <- argv[4]

path <- "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/exp" 
design_file <- "input_edgeR_ColDesign.txt"
data_matrix_file <- "brca_rnaseq_rawCount_l3_geneSymbol.mat"
output_file <- "brca_rnaseq_rawCount_l3_geneSymbol" 
outSampleFile <- "final.header_exp.sorted"

##----------------------------
#load data
dataMat <- read.delim(data_matrix_file,header=T)
header <- as.character(unlist(read.delim(data_matrix_file,header=F,nrow=1)[-1]))
# dataMat[grep("SLC35E2",dataMat[,1])[2],1] = as.factor("SLC35E2B")
rnames  = as.character(dataMat[,1])
rnames[16273] = "SLC35E2B"
rownames(dataMat) <- rnames
dataMat <- dataMat[,-1]

#load desgin file
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

# ##----------------------------
# #fitting model
library(edgeR)
setwd(path)
y <- DGEList(counts=dataMat,group=condition)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
save.image(file=paste(output_file,".rda",sep=""))


###########----------------------------
#output --- plot
pdf(paste(output_file,".pdf",sep=""))
plotBCV(y, cex=0.4)

de <- decideTestsDGE(et, p=0.05, adjust="BH")
detags <- rownames(topTags(et, n=500))
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "blue")

outTag 		<- topTags(et,n=Inf,adjust.method="BH",sort.by="p.value")
topn 		= length(rownames(outTag)[outTag@.Data[[1]]$FDR < 0.01])
# rownames(outTag)[outTag@.Data[[1]]$FDR < 0.01]
outTag 		<- topTags(et,n=topn,adjust.method="BH",sort.by="p.value")

write.table(outTag,file=paste(output_file,"_DEG_stats.txt",sep=""),sep="\t",quote=F)
require(DESeq)
require(vsn)
cds    		<- newCountDataSet(dataMat[rownames(outTag),],condition)
cds   		<- estimateSizeFactors( cds )
cdsBlind 	= estimateDispersions(cds, method = "blind")
vsd 		<- getVarianceStabilizedData(cdsBlind)
outTagMat 	<- getVarianceStabilizedData(cdsBlind)

par(mfrow=c(2,1))
notAllZero = (rowSums(counts(cds))>0)
meanSdPlot(log2(counts(cds)[notAllZero, ] + 1), ylim = c(0,2.5))
meanSdPlot(outTagMat[notAllZero, ], ylim = c(0,2.5))

dev.off()

###output file------------
colnames(outTagMat) 	<- header
outSample 				<- as.character(unlist(read.table(outSampleFile,header = F)))
outSampleIdx 			<- sapply(outSample,function(x){grep(x,colnames(outTagMat))})
out 					<- outTagMat[,outSampleIdx]
outTagMat_header		<- as.matrix(t(c("Gene",colnames(out))))

write.table(outTagMat_header,file=paste(output_file,"_DEG_Mat.txt",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
write.table(out,file=paste(output_file,"_DEG_Mat.txt",sep=""),sep="\t",quote=F,append = TRUE,col.names=F)
