library(affy)
setwd("/ifs/data/c2b2/ac_lab/jh3283/AML/Valk_AML/CEL/")
expdata <- ReadAffy()

eset <- rma(expdata)
temp <- eset
eset2 <-new("ExpressionSet", exprs = as.matrix(temp), annotation="hgu133aALIAS2PROBE")
write.exprs(eset, file="/ifs/data/c2b2/ac_lab/jh3283/AML/Valk_AML/valk_aml_affyu133_293.exp")
mypm <- pm(expdata)
mymm >- mm(expdata)




library(affyQCReport)
QCReport(expdata, file="valk_aml_QC.pdf") # Generates a comprehensive QC report for the AffyBatch object 'mydata' in PDF format. See affyQCReport for details.
# deg <- AffyRNAdeg(expdata); 
# summaryAffyRNAdeg(deg); 
# plotAffyRNAdeg(deg) # Performs RNA degradation analysis. It averages on each chip the probes relative to the 5'/3' position on the target genes. A summary list and a plot are returned.
# image(mydata[ ,1]) # Reconstructs image with log intensities of first chip.
# hist(mydata[ ,1:2]) # Plots histogram of PM intensities for 1st and 2nd array.
# hist(log2(pm(mydata[,1])), breaks=100, col="blue") # Plos bar histogram of the PM ('pm') or MM ('mm') log intensities of 1st array.
par(mfrow=c(1,1))
pdf("qc_boxplot_inten.pdf")
boxplot(expdata,col="red") # Generates a box plot of un-normalized log intensity values.
boxplot(data.frame(exprs(eset)),col="blue", main="Normalized Data") # Generates a box plot of normalized log intensity values.
dev.off()
mva.pairs(pm(mydata)[,c(1,4)]) # Creates MA-plot for un-normalized data. A MA-plot is a plot of log-intensity ratios (M-values) versus log-intensity averages (A-values) between selected chips (here '[1,4]').
mva.pairs(exprs(eset)[,c(1,4)]) # Creates MA-plot for normalized data.


expMA <- exprs(eset)
colnames(expMA) <- gsub(".CEL.gz","",colnames(expMA))

samples <- colnames(expMA)
NDsample <- c("GSM20684","GSM20714","GSM20764","GSM20767","GSM20772","GSM20775","GSM20777","GSM20778","GSM20781","GSM20797","GSM20816","GSM20820","GSM20821","GSM20829","GSM20848","GSM20885","GSM20887","GSM20888","GSM20892","GSM20899","GSM20900","GSM20952","GSM20962")

exp270MA <- expMA[,setdiff(samples,NDsample)]
design270<- read.table("/ifs/home/c2b2/ac_lab/jh3283/schome/AML/valk/valk_aml_design_270.txt",sep = "\t",header=F)
row.names(design270) <- design270[,1]
design270 <- design270[,2:271]
write.table(exp270MA,"valk_aml_270.exp")

sM012 <- design270$ID_REF[which(design270$FAB_stage==1)]
sM3 <- design270$ID_REF[which(design270$FAB_stage==2)]
sM456 <- design270$ID_REF[which(design270$FAB_stage==3)]
scontrol<- design270$ID_REF[which(design270$FAB_stage==0)]



myttest <- function(exp.matrix,group1,group2) {
	rownum <- length(exp.matrix[,1])
	tprobDF <- data.frame(probe=row.names(exp.matrix),t=rep(0,rownum),pvalue=rep(0,rownum))	

	for (i in 1: rownum){
		tprobDF$t[i] <-t.test(exp.matrix[i,group1],exp.matrix[i,group2])$statistic
		tprobDF$pvalue[i] <-t.test(exp.matrix[i,group1],exp.matrix[i,group2])$p.value
	}
	return(tprobDF)
}




getGenelist <- function(tresult,g_num=0,p_cutoff=0){
	glist <- list()
	if(p_cutoff > 0)glist <- tresult$probe[which(tresult$pvalue < p_cutoff)]
    if(g_num > 0) glist <- tresult$probe[order(tresult$pvalue)][1:g_num] 
    return(glist)
}




### annotation
library("hgu133a.db")

# myAnnot <- data.frame(ACCNUM=sapply(contents(hgu133aACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(hgu133aALIAS2PROBE), paste, collapse=", "), DESC=sapply(contents(hgu133aGENENAME), paste, collapse=", "))

# contents(hgu133aALIAS2PROBE)[1:20]

myAnnot <- function(geneList){
	len <- length(geneList)
	glist.new <- data.frame(probe=geneList,SYMBOL= rep(NA,len))
	for (i in 1:len){
		tmpg <- as.character(glist.new$probe[i])
		glist.new$SYMBOL[i] <- get(tmpg,env=hgu133aSYMBOL)
	}
	return(glist.new)
}


 tresultM012 <- myttest(exp270MA,sM012,scontrol)
 tresultM3 <- myttest(exp270MA,sM3,scontrol)
 tresultM456 <- myttest(exp270MA,sM456,scontrol)

 geneLM012 <- getGenelist(tresultM012,p_cutoff=0.0001)
 geneLM3 <- getGenelist(tresultM3,p_cutoff=0.0001)
 geneLM456 <- getGenelist(tresultM456,p_cutoff=0.0001)

  geneLM012 <- getGenelist(tresultM012,g_num=250)
 geneLM3 <- getGenelist(tresultM3,g_num=250)
 geneLM456 <- getGenelist(tresultM456,g_num=250)


geneSybmolM3 <- myAnnot(geneLM3)
geneSybmolM456 <- myAnnot(geneLM456)
geneSybmolM012 <- myAnnot(geneLM012)




write.table(na.omit(geneSybmolM012$SYMBOL),"valk_270_GL_m012.txt",col.names=F,row.names=F)
write.table(na.omit(geneSybmolM3$SYMBOL),"valk_270_GL_m3.txt",col.names=F,row.names=F)
write.table(na.omit(geneSybmolM456$SYMBOL),"valk_270_GL_m456.txt",col.names=F,row.names=F)



geneSybmolM3 <- myAnnot(geneLM3)
geneSybmolM456 <- myAnnot(geneLM456)
geneSybmolM012 <- myAnnot(geneLM012)




write.table(na.omit(geneSybmolM012$SYMBOL),"valk_270_GL_m012.txt",col.names=F,row.names=F)
write.table(na.omit(geneSybmolM3$SYMBOL),"valk_270_GL_m3.txt",col.names=F,row.names=F)
write.table(na.omit(geneSybmolM456$SYMBOL),"valk_270_GL_m456.txt",col.names=F,row.names=F)




# get("200657_at" ,hgu133aSYMBOL)
# get("200657_at" ,hgu133aACCNUM)




