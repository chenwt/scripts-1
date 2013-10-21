#wd /ifs/scratch/c2b2/rr_lab/akl2140/projCompMethods/wd/feature
#!/bin/Rscirpts
#Author: Jing He
#Date: 
#Last Update: 
#parameters: 
#example: Rscript col2BinMatrix.r ICD9

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
data <- read.table(args,header=T)
print(paste("input data frame of", dim(data)[1],"*",dim(data)[2]))
colnames(data) <- c("id","var")

# out <- model.matrix(~var,data=data,contrasts.arg=list(var=diag(nlevels(data$var))))[,-1]
out <- sparse.model.matrix(~var,data=data,contrasts.arg=list(var=diag(nlevels(data$var))))[,-1]

out <- cbind(data$id,as.data.frame(out))
colnames(out) <- c("id",levels(data$var))
out <- ddply(out,.(id),numcolwise(sum))

# data <- apply(data,2,function(x){as.factor(as.character(x))})
# out <- modle.matrix(~ .,data=data[,-1],
# 		contrasts.arg=lapply(data[,-1],contrasts,contrasts=FALSE))
# out <- ddply(out,.(id),numcolwise(sum))
write.table(out,paste(args,".matrix",sep=""),sep="\t",quote=FALSE,row.names=T,col.names=T)
print(paste("output matrix of", dim(out)[1],"*",dim(out)[2],"in", paste(args,".mtx",sep="")))
# for (i in 1:ceiling(dim(outMX)[1]/30000){
# 	write.matrix(outMX[((i-1)*30000 + 1):min(i*30000,dim(outMX)[1])],"FEATURE-ICD9.simplified.matrix",sep="\t")
# }

