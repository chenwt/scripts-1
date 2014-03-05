##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep.26,13 
#Last Updated:
#Usage:
#Description:plot script for labmeeting presentation
#source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projAML/WXS/plot_MrMtgPreppi.r")

setwd("/Volumes/ac_lab/jh3283/SCRATCH/projAML/WXS/reports/result_labmetgSep24") 

dstat <- read.table("MR_mutgene_preppi_stat.txt",header=F)
dcorr <- read.table("MR_mutgene_preppi.txt",sep="\t")
df <- dstat[which(dstat$V6>0),]
df <- df[order(df$V1),]


require(reshape2)
# df$V6 <- as.character(df$V6)
df_wide <- dcast(
	transform(df,id=rep(1:length(unique(df$V1)),times=table(df$V1))),
	id ~ V3, value.var="V6",fill=0)
df_wide <- df_wide[,-1]
df_wide <- apply(df_wide,2,as.numeric)
row.names(df_wide) <- names(table(df$V1))
pdf("fig/fig_heatmap_mr_mutgene_stat.pdf")
heatmap(df_wide,labRow=row.names(df_wide),
	col=colorRampPalette(c("white","blue"))(256),
	xlab="mutated genes",Rowv=NA,Colv=NA)
dev.off()


# m <- df_wide
# colnames(m) <- colnames(m)
# rownames(m) <- row.names(m)
# # m

# image(1:ncol(m), 1:nrow(m), t(m), col = colorRampPalette(c("white","blue"))(256), axes = FALSE)
# axis(1, 1:ncol(m), colnames(m))
# axis(2, 1:nrow(m), rownames(m))
# for (x in 1:ncol(m))
#   for (y in 1:nrow(m))
# 	  if(m(y,x) >0 )
# 	    text(x, y, m[y,x])



