##!/usr/bin/Rscript
#Author: Jing He
#Date:Oct 3,2013 
#Last Updated:
#Usage:Rscript do_diffExprDESeq.r <input1> <input2> <input3>
#input: 	<read fpkm matrix> 
# 			<BED file which generate the matrix above, row annotation>
# 			<sample design file, column annotation>
# out: <diff expreseed lincRNAs and their expression>
# source("~/scripts/projRe/do_diffExprFPKM.r")
#  get big FPKM matrix getFPKM.sh


wd <- "/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/"
setwd(wd)
#read data

inf1 <- "result_genes_fpkm_matrix_17"
# inf2 <- "/ifs/scratch/c2b2/ac_lab/jh3283/database/ncRNAlib/broad/lincRNAs_transcripts.bed"
inf2 <- "/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/FPKM.trackid.annot"
inf3 <- "/ifs/scratch/c2b2/ac_lab/jh3283/projRe/brca/result/barcode.list.fpkm"
fpkmdf <- read.table(inf1)
rns <- fpkmdf[,1]; fpkmdf <- fpkmdf[,-1]
rowannot <- read.table(inf2,sep="\t")
colannot <- read.table(inf3)
rownames(fpkmdf) <- fpkmdf$V1
colnames(fpkmdf) <- colannot$V1

#extract different condition

get_grps <- function(colannot){
	colgrps <- sapply(unique(colannot$V2),function(x)which(colannot$V2 == x))
	names(colgrps) <- unique(colannot$V2)
	return(colgrps)
}

fpkmData <- fpkmdf
do_Expr = function(cond1,cond2,fpkmData,colannot){
	colgrps <- get_grps(colannot)
	idx <- c(unlist(colgrps[[cond1]]),unlist(colgrps[[cond2]]))
	# print(unique(colannot))
	condition <- factor(colannot$V2[idx])
	fpkmTuNo <- fpkmData[,idx]
	rownames(fpkmTuNo) <- rownames(fpkmData)
	#QC
	print("QC...")
	fpkmd <- 
	cds = estimateSizeFactors( cds )
	# fpkms( cds, normalized=TRUE )
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


# cond1 <- "TM" 
cond1 <- "TP" 
cond2 <- "NT"
# res <- do_Expr(cond1,cond2,fpkmdf,colannot)
# res$chr_pos = rns
# res$annot = rowannot$V4
# outfpkm = 1
# write.csv(res, file=paste("result",cond1,"_",cond2,"_",outfpkm,".csv",sep=""),row.names=F,quote=F )