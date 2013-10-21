#!/bin/Rscirpts
#Author: Jing He
#Date: Mar.24,2013
#Last Update: Mar.24,2013

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
  #   print(snp$Altfreq[((i-1)*bin+1):(i* bin)])
  #   print (" ") 
  # }
  modes[modes<0.2] <- 0
  modes[modes>0.8] <- 1
  ff <- unlist(strsplit(freq.file,"/"))
  tt <- gsub(".freq","",ff[length(ff)])
  plot(snp$Pos[1:(bin-1)],modes[,1],ylim=c(0,1),
              type="l",col="blue",main=tt, xlab="",ylab="",xaxt="n",yaxt="n")
  lines(snp$Pos[1:(bin-1)],modes[,2],ylim=c(0,1),
              type="l",col="blue",xlab="",ylab="",xaxt="n",yaxt="n")
}


plotCNV22<- function(freq_mode.file, snp.file, out.name) {
###------------------------------ load individual snp
  require(ggplot2)
  require(grid)
  require(scales)
  
  pdf(out.name,width=18, height=12)
  # grid.newpage()
  # pushViewport(viewport(layout=grid.layout(6, 4)))
  par(mfrow=c(4,6),mai=c(0.1,0.3,0.2,0.1))
  for (i in 1:22) {
    # row <- floor((i - 1) / 4) + 1
    # column <- (1 + (i + 3) %% 4)
    freq.file <- gsub("CHR", i, freq_mode.file)
    print(getCNV(freq.file,snp.file)
      # , vp=vplayout(row, column)
      )
  }
  # grid.text("test", vp=viewport(layout.pos.row = 6, layout.pos.col=3:4))
  dev.off()

}


#get the two mode of a vector
myMode <- function(v){
  # require(KernSmooth)
  # require(diptest)
  v <- as.numeric(v)
   require(modeest)
   return(c(m1=round(mlv(v[v<0.5],method="lientz",bw=0.01)$M,digits=2), 
      m2=round(mlv(v[v>=0.5],method="lientz",bw=0.01)$M,digits=2)))

  # print(plot(bkde(v),col="blue"),type="l",xlab="",ylab="")
  # tp <- dip.test(bkde(v)$y)$p.value
  # br <- seq(0, 1, by=0.1)
  # vt <- summary(cut(v, br, labels=br[-2], include.lowest=T, ordered_result=T))
  # if(tp < 0.01){ # test for bimodal distribution
  #   m1 <- min(as.numeric(names(vt[order(vt,decreasing=T)])[1:2])) 
  #   m2 <- max(as.numeric(names(vt[order(vt,decreasing=T)])[1:2])) 
  #   # print(c(m1,m2)) 
  #   return(c(m1,m2))
  # }else if(tp >= 0.01)
  # {
  #   m <- as.numeric(names(vt[order(vt,decreasing=T)])[1])
  #   return(c(m,m))}



}