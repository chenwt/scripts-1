##!/usr/bin/Rscript
#Author: Jing He
#Date:Oct 3,2013 
#Last Updated:
#Usage:Rscript do_diffExprDESeq.r <input1> <input2> <input3>
#input: 	<read count matrix> 
# 			<BED file which generate the matrix above, row annotation>
# 			<sample design file, column annotation>
# out: <diff expreseed lincRNAs and their expression>
# source("~/scripts/projRe/do_diffExprDESeq.r")

wd <- "/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/"
setwd(wd)
#read data
inf1 <- "result_tss_count_17.txt"
inf2 <- "/ifs/scratch/c2b2/ac_lab/jh3283/database/ncRNAlib/broad/lincRNAs_transcripts.bed"
inf3 <- "barcode.list.sample_code"
countdf <- read.table(inf1)
rns <- countdf[,1]; countdf <- countdf[,-1]
rowannot <- read.table(inf2,sep="\t")
colannot <- read.table(inf3)
# rownames(countdf) <- countdf$V1
colnames(countdf) <- colannot$V1

#extract different condition

get_grps <- function(colannot){
	colgrps <- sapply(unique(colannot$V2),function(x)which(colannot$V2 == x))
	names(colgrps) <- unique(colannot$V2)
	return(colgrps)
}

do_model = function(cond1,cond2,countData,colannot){
	require(DESeq)
	colgrps <- get_grps(colannot)
	idx <- c(unlist(colgrps[[cond1]]),unlist(colgrps[[cond2]]))
	# print(unique(colannot))
	condition <- factor(colannot$V2[idx])
	countTuNo <- countData[,idx]
	# rownames(countTuNo) <- rownames(countData)
	#QC
	print("QC...")
	cds <- newCountDataSet(countTuNo,condition)
	cds = estimateSizeFactors( cds )
	# counts( cds, normalized=TRUE )
	cds = estimateDispersions( cds )
	# str( fitInfo(cds) )
	figcount = 1
	pdf(paste("plot_DESeq_",cond1,"_",cond2,"_",figcount,".pdf",sep=""))
	plotDispEsts(cds)
	dev.off()

	#diff exp
	print("diff...")
	res = nbinomTest( cds, cond1, cond2 )
	figcount = figcount + 1
	pdf(paste("plot_DESeq_",cond1,"_",cond2,"_",figcount,".pdf",sep=""))
	plotMA(res)
	hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
	dev.off()
	return(res)
	print("Done!")
}


cond1 <- "TM" 
# cond1 <- "TP" 
cond2 <- "NT"
colgrps <- get_grps(colannot)
res <- do_model(cond1,cond2,countdf,colannot)
res$chr_pos = rns
res$annot = rowannot$V4
outcount = 1
write.csv(res, file=paste("result",cond1,"_",cond2,"_",outcount,".csv",sep=""),row.names=F,quote=F )


# resSig = res[ res$padj < 0.1, ]
# resSig$annot = rowannot[which(res$padj < 0.1,ind.)]

