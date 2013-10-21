#!/bin/Rscirpts
#Author: Jing He
#Date: 
#Last Update: 
#parameters: 
#example: Rscript col2SparseMtx.r ICD9.simplified
#wd /ifs/scratch/c2b2/rr_lab/akl2140/projCompMethods/wd/feature

###------------------------------header------------------------------
args <- commandArgs(TRUE)
if (is.null(args)){
	print("Please provide parameters")
	exit
}else{
	print(paste("Input file:",args,sep=""))
}
###------------------------------start coding##------------------------------
require("plyr")
require("Matrix")
require("MASS")
require("glmnet")
require("reshape2")

dataraw <- read.table(args,header=T)
# dataraw <- read.table("FEATURE-ICD9.simplified.100patients.ICD9")
colnames(dataraw) <- c("pid","var")
#dataraw[order(dataraw$pid),]
dataraw$value <- 1

print(paste("input dataraw frame of", dim(dataraw)[1],"*",dim(dataraw)[2]))

dataC <- dcast(dataraw,pid ~ var + value, length, value.var="value")

dataMX <- Matrix(as.matrix(dataC),sparse=TRUE)

print("print image")
#pdf(paste(args,".pdf",sep=""))
image(dataMX, main = "Patients ICD9")
dev.off()
# col.names(out) <- col.names(data)
# out <- sparse.model.matrix(~ .,data=data,contrasts.arg=lapply(data[,-1],contrasts,contrasts=FALSE), row.names=T)[,-1]
# write.table(out,paste(args,".matrix",sep=""),sep="\t",quote=FALSE,row.names=T,col.names=T)

print("writing out the data")
save(dataMX,file=paste(args,".Rdata",sep=""))
print("finished")
print(paste("save sparse matrix of", dim(out)[1],"*",dim(out)[2],"in", paste(args,".Rdata",sep="")))


# image(out,
# xlim = .5 + c(0, di[2]),
# ylim = .5 + c(di[1], 0), aspect = "iso",
# sub = sprintf("Dimensions: %d x %d", di[1], di[2]),
# xlab = "Column", ylab = "Row", cuts = 15
# )
