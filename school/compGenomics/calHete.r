#!/bin/Rscirpts
#Author: Jing He
#Date: Mar, 2013
#Last Update: 

# calculate the heterogeneity level for each snp using tumor v.s. control 
getHeteCHR <- function(t.file, n.data, snp.file, bin.num){
	tsnp <- getSNP(t.file, snp.file)
	nsnp <- getSNP(n.file, snp.file)
	mtn <- merge(nsnp,tsnp, by=c("Chr","Pos")) 
	# the filtered out snps might be interesting because they only appeared in tumor sample
	apply(mtn[1:5,c(9,5)],1,mySNPHete)

}

## get the deletion region using chromosome freq data
getCNV <- function(freq.file,snp.file) {
	snp <- getSNP(freq.file, snp.file)
	snp <- snp[order(snp$Pos),]
	len <- nrow(snp)
	bin <- round(len / 2)
	from <- 1
	to <- bin
	# hist(snp$Altfreq[from:to])
	modes <- t( vapply(1:(bin-1), FUN=function(x)
				{myMode(as.numeric(snp$Altfreq[x:(x+bin-1)]))} ,
				vector(mode="numeric",length=2)) )

	# for (i in 1:bin.num){
	# 	print(snp$Altfreq[((i-1)*bin+1):(i* bin)])
	# 	print (" ") 
	# }
	modes[modes<0.2] <- 0
	modes[modes>0.8] <- 1
	ff <- unlist(strsplit(freq.file,"/"))
	tt <- gsub(".freq","",ff[length(ff)])
	print(plot(snp$Pos[1:(bin-1)],modes[,1],type="l",ylim=c(0,1),col="blue",main=tt))
    lines(snp$Pos[1:(bin-1)],modes[,2],type="l",ylim=c(0,1),col="blue")
}

plotCNV22<- function(freq_mode.file, snp.file, out.name) {
###------------------------------ load individual snp
	require(ggplot2)
	require(grid)
	require(scales)
  
  pdf(out.name)
  # grid.newpage()
  # pushViewport(viewport(layout=grid.layout(6, 4)))
  par(mfrow=c(6,4))
  for (i in 1:22) {
    row <- floor((i - 1) / 4) + 1
    column <- (1 + (i + 3) %% 4)
    freq.file <- gsub("CHR", i, freq_mode.file)
    print(getCNV(freq.file,snp.file)
    	# , vp=vplayout(row, column)
    	)
  }
  # grid.text("test", vp=viewport(layout.pos.row = 6, layout.pos.col=3:4))
  dev.off()

}
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)






###------------------------------ get SNP data
getSNP <- function(freq.file,snp.file){
	freq <- read.table(freq.file, skip=5, header=TRUE)
	freq$Altfreq <- freq$Altreads / freq$Totalreads
	freq <- freq[freq$Altfreq > -1 & freq$Altfreq < 2, ]
	snp <- read.table(snp.file,sep="\t",header=F)
	colnames(snp) <- c("Chr","Pos","Ale")
	#------------------------------ load 22 freq data
	snpFreq <- merge(x=freq,y=snp,by=c("Chr","Pos"))
	snpFreq$Altfreq <- snpFreq$Altreads / snpFreq$Totalreads
	return(snpFreq)
}



myBinomTest <- function(x,n,p,alpha = 0.05){
 return(binom.test(x,n,p,alternative="two.sided")$p.value < alpha)
}

#get the two mode of a vector
myMode <- function(v){
	# require(KernSmooth)
	# require(diptest)
	v <- as.numeric(v)
	# print(plot(bkde(v),col="blue"),type="l",xlab="",ylab="")
	# tp <- dip.test(bkde(v)$y)$p.value
	# br <- seq(0, 1, by=0.1)
	# vt <- summary(cut(v, br, labels=br[-2], include.lowest=T, ordered_result=T))
	# if(tp < 0.01){ # test for bimodal distribution
	# 	m1 <- min(as.numeric(names(vt[order(vt,decreasing=T)])[1:2])) 
	# 	m2 <- max(as.numeric(names(vt[order(vt,decreasing=T)])[1:2])) 
	# 	# print(c(m1,m2)) 
	# 	return(c(m1,m2))
	# }else if(tp >= 0.01)
	# {
	# 	m <- as.numeric(names(vt[order(vt,decreasing=T)])[1])
	# 	return(c(m,m))}

	 require(modeest)
	 return(c(m1=round(mlv(v[v<0.5],method="lientz",bw=0.01)$M,digits=2), 
	 		m2=round(mlv(v[v>=0.5],method="lientz",bw=0.01)$M,digits=2)))


}

# given the BAF of a SNP, calculate the heterogeneity level
mySNPHete <- function(x){
	t <- x[1];n <- x[2]
	if(n < 0.5 & n > 0 & t > 0){
		return(2- 1/t)}
	else if(n >= 0.5){
		return((1-2 * t)/(1-t))
	}
}


PlotHete <- function(freq.file,snp.file,bin.num) {
	SNPdata <- getSNP(freq.file,snp.file)
	SNPdata <- SNPdata[order(SNPdata$Pos),]
	# bin.size <- (as.numeric(max(SNPdata$Pos)) - as.numeric(min(SNPdata$Pos)) / as.numeric(bin.num))
	from <- 1
	to <- bin.size
	# hist(SNPdata$Altfreq[from:to])
	bin.size <- 10000000
	ppos <- seq(min(SNPdata$Pos), max(SNPdata$Pos), by=bin.size)
	row.names(SNPdata) <- SNPdata$Pos
	mode <- vapply(c(1,ppos), FUN=function(x){myMode(unlist(SNPdata$Altfreq[x:(x+bin.size-1)]))},vector(mode="numeric",length=2))
	myMode(SNPdata$Altfreq[1:ppos[1]])
	plot(y=mode2[1],x=ppos)
	lines(y=mode2[2],x=ppos)
	# x <- as.numeric(SNPdata$Altfreq)[from:to]
	# br <- seq(0, 1, by=0.1)
	# xt <- summary(cut(x[x!=0], br, labels=br[-1], include.lowest=T, ordered_result=T))
	# xt[order(xt,decreasing=T)]
}
