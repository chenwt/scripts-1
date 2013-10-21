##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep, 20,2013 
#Last Updated:
#Usage:
#Description: do analysis to the table generated 

rootd <- "~/SCRATCH/projAML/WXS/reports"
setwd(rootd)
data <- unique(read.table("result_varFreq_v2_overlap.txt",header=T,fill=NA))

genes <- unique(data$gene_name)
data$nat <- as.integer(data$N_Alt) / as.integer(data$N_Total)
data$tat <- as.integer(data$T_Alt) / as.integer(data$T_Total)
data$rat <- as.integer(data$R_Alt) / as.integer(data$R_Total)
x = genes[1]

plotGeneMAF <- function(x){
	plotdf <- data[which(data$gene_name==x),c("nat","tat","rat")]
	rownames(plotdf) <- paste(data[which(data$gene_name==x),"chr_pos"],
							data[which(data$gene_name==x),"PID"],sep="_")
	png(paste("graph/",x,".png",sep=""))
	plot(0,0,xlim = c(0,4),ylim = c(min(plotdf),min(max(plotdf),1)),
		type = "n",xlab=x，ylab="MAF")
	apply(plotdf,1,function(row){lines(x=1:3,y=row,type="b")})
	dev.off()
}

sapply(genes,plotGeneMAF)