#!/bin/Rscirpts
#Author: Jing He
#Date: Mar. 2013
#Last Update: Mar.23,2013

getdata <- function(filename){
  Freqchr <- read.table(filename,skip=5,header=T,sep="\t")
  Freqchr <- Freqchr[Freqchr$Altreads > -1 & Freqchr$Altreads/Freqchr$Totalreads < 2,]
  Freqchr$Altfreq <- Freqchr$Altreads/Freqchr$Totalreads
  return(Freqchr)
}

###------------------------------
#get snp freqency of one chromosome according to exsiting snp file and freq file
getSnpFreq <- function(filename,snp.file){
  snp <- read.table(snp.file,sep="\t",header=F)
  colnames(snp) <- c("Chr","Pos","Ale")
  snp <- snp[-grep("X|Y",snp$Chr),]
  snp$Chr <- as.integer(snp$Chr)
  Freq <- getdata(filename)
  Freqchr <- merge(x=snp,y=Freq,by=c("Chr","Pos"))
  Freqchr$Altfreq <- Freqchr$Altreads / Freqchr$Totalreads
  return(Freqchr)
}

getSnpFreq22Data <- function(filename,snp.file){
  # source("/Volumes/ys_lab/jh3283/scripts/net/PlotFreq.r")
  ###------------------------------ load individual snp
  snp <- read.table(snp.file,sep="\t",header=F)
  colnames(snp) <- c("Chr","Pos","Ale")
  snp <- snp[-grep("X|Y",snp$Chr),]
  snp$Chr <- as.integer(snp$Chr)

  ###------------------------------ load 22 freq data
  Freq <- data.frame(Chr=integer(),Pos=integer(),Totalreads=integer(),Altreads=integer())
  for (i in 1:22){
    filename <- filename
    filename <- gsub("CHR",i,filename)
    Freq <- rbind(Freq,getdata(filename))
  }
  Freqchr <- merge(x=snp,y=Freq,by=c("Chr","Pos"))
  Freqchr$Altfreq <- Freqchr$Altreads / Freqchr$Totalreads
  # Freqchr$Zscore <- (Freqchr$Altfreq - mean(Freqchr$Altfreq))/sd(Freqchr$Altfreq)
  return(Freqchr)
}

t.TestHetegeneity <- function(ac.num, normal_mode.file, snp.file) {

    tumor <- getSnpFreq22Data(freq_mode.tumor.file,snp.file)
    normal <- getSnpFreq22Data(freq_mode.normal.file, snp.file)



} 

# myBinomialTest <- function(){
  
# }

###------------------------------ 
## case-control test(x-y), for each snp
# input x/y.vector AltAllele freq,total allele freq
# 

myFETest <- function(x,y){
 x[2] <- x[2] - x[1]
 y[2] <- y[2] - y[1]
 z <- rbind(x,y)
 res.t <- fisher.test(z)
 res<- vector(mode="numeric",length=2)
 res <- c(res.t$estimate,res.t$p.value)
 names(res) <- c("OR.es","p.value")
 # if(res[1] < 1){
 #  res[1] <- (1/res[1])
 # }
 return(res)
}


###------------------------------
# two freq file after call BAF
# fixed normal file 
# snp.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt"

myChrFETest <- function(ac.num, chr.num, snp.file){
  ###------------------------------load scripts and file
  source("~/Dropbox/scripts/projNET/BAFTest.r")
  freq_mode.tumor.file <- "/Volumes/ys_lab_scratch/jh3283/net/ACACCODE/BAF/chrCHR.freq"
  freq_mode.tumor.file <- gsub("ACCODE",ac.num,freq_mode.tumor.file)
  freq_mode.normal.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chrCHR.freq"

  freq.tumor.file <- gsub("CHR",chr.num,freq_mode.tumor.file)
  freq.normal.file <- gsub("CHR",chr.num,freq_mode.normal.file)

  t <- getdata(freq.tumor.file)
  n <- getdata(freq.normal.file)
  # print("file loaded\n")
  ### fisher exact test
  # print("start merging...\n")
  freq.common <- merge(t,n,by=c("Chr","Pos"),suffixes=c("T","N"))
  # print("start testing....\n")
  res.test <- t(apply(freq.common,1,function(x){myFETest(x[c(4,3)],x[c(7,6)])}))
  rownames(res.test) <-freq.common$Pos
  ## p value adjustment Holm–Bonferroni method
  res.test <- data.frame(res.test)
  res.test <- res.test[order(res.test$p.value),]
  alpha <- 0.05
  M <- nrow(res.test)
  res.test <- cbind(res.test,alpha.adj = vapply(1:length(res.test$p.value),FUN=function(x){alpha/(M+1-x)},1))
  res.test <- res.test[apply(res.test,1,function(x){
                                          if(x[2] > x[3])return(FALSE)
                                          else return(TRUE)
                                          }),]

  return(res.test)
}

#------------------------------ to all chromosome
my22ChrFETest <-function(ac.num,snp.file) {
  source("~/Dropbox/scripts/projNET/BAFTest.r")
  meta <- data.frame(OR.es=numeric(),p.value=numeric(),alpha.adj=numeric())
  res <- data.frame(Chr=numeric(),Pos=numeric(),OR.es=numeric(),p.value=numeric(),alpha.adj=numeric())

  for (i in 1:22){
    print(paste("analysising chromosome",i))
    meta <- myChrFETest(ac.num,i,snp.file)
    meta$Chr <- i
    meta$Pos <- rownames(meta)
    meta <- meta[,c("Chr","Pos","OR.es","p.value","alpha.adj")]
    res <- rbind(res,meta)
  }
  return(res)
}


MAplot <- function(tumor.file,normal.file,t.num,n.num,del_reg.file){
  require(ggplot2)
  ###------------------------------ load data
  tumor <- read.table(tumor.file,sep="\t",header=T,stringsAsFactors=F)
  normal <- read.table(normal.file,sep="\t",header=T,stringsAsFactors=F)
  tumor$FPKM[tumor$FPKM == 0] <- 1e-25
  normal$FPKM[normal$FPKM == 0] <- 1e-25
  M <- log(tumor$FPKM,2)-log(normal$FPKM,2)
  A <- 1/2 * (log(tumor$FPKM,2) + log(normal$FPKM,2))
  data.plot <- data.frame(M=M,A=A)
  M[M > 0]
  ###------------------------------load deletion region data
  del.reg <- read.table(del_reg.file,sep="\t",header=T)
  del.reg <- del.reg[,2:4]
  colnames(del.reg) <- c("chromosome","start","end")
  del.reg$chromosome <- as.numeric(del.reg$chromosome)
  tmpvar <-  unlist(lapply(tumor$locus,FUN = function(x){unlist(strsplit(x,"-"))[1]}))
  emph.chr <-  unlist(lapply(tmpvar,FUN = function(x){unlist(strsplit(x,":"))[1]})) 
  emph.start <- as.numeric( unlist(lapply(tmpvar,FUN = function(x){unlist(strsplit(x,":"))[2]})) )
  emph.end <- as.numeric( unlist(lapply(tumor$locus,FUN = function(x){unlist(strsplit(x,"-"))[2]})))
  emph <- cbind(data.plot, emph.chr,emph.start,emph.end)
  colnames(emph) <- c(colnames(data.plot),"chr","start","end")

  ##------------------------------
  emph.plot <- data.frame(M=numeric(),A=numeric(),chr=numeric(),start=numeric(),end=numeric())
  emph$chr <- as.numeric(emph$chr)

  for (i in 1:nrow(emph)) {
    for (j in 1:nrow(del.reg)){
      if( !is.na(emph[i,3]) & emph[i,3] == del.reg[j,1]){
        if( emph[i,4]>= del.reg[j,2] & emph[i,5] <= del.reg[j,3]){
              emph.plot <- rbind(emph.plot,emph[i,])
            break
        }
      }
    next
    }
  }


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
  
  pdf(out.name,width=15, height=12)
  # grid.newpage()
  # pushViewport(viewport(layout=grid.layout(6, 4)))
  par(mfrow=c(4,6),mai=c(0.1,0.2,0.1,0.1))
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
  #   m1 <- min(as.numeric(names(vt[order(vt,decreasing=T)])[1:2])) 
  #   m2 <- max(as.numeric(names(vt[order(vt,decreasing=T)])[1:2])) 
  #   # print(c(m1,m2)) 
  #   return(c(m1,m2))
  # }else if(tp >= 0.01)
  # {
  #   m <- as.numeric(names(vt[order(vt,decreasing=T)])[1])
  #   return(c(m,m))}

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



  
