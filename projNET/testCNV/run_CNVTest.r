#!/bin/Rscirpts
#Author: Jing He
#Date: 
#Last Update: 
#parameters: 
#example: Rscripts ~/scripts/projNET/run_BAFTest.r 13 19 AC1SNPs.txt

###------------------------------header------------------------------
args <- commandArgs(TRUE)

if (is.null(args)){
  print("Please provide parameters")
  exit
}else{
  print(args)
}

###------------------------------start coding##------------------------------

args <- commandArgs(TRUE)

t.num <- 13  # 10,13,16 
n.num <- 19
snp.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt"
Sys.setenv( NETDIR="/Volumes/ys_lab_scratch/jh3283/net/",
      SCRIPTS="/Users/jh3283/Dropbox/scripts/"
      )
NETDIR <- Sys.getenv("NETDIR")
CWD <- paste(NETDIR,"AC",t.num,"/res/",sep="")
setwd(CWD)
source("~/Dropbox/scripts/projNET/testCNV/CNVTest.r")
# freq.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC13/BAF/chr1.freq"
freq_mode.tumor <- paste(NETDIR,"AC",t.num,"/BAF/","chrCHR.freq",sep="")
freq_mode.normal <- paste(NETDIR,"AC",n.num,"/BAF/","chrCHR.freq",sep="")

freq.tumor <- gsub("CHR",1,freq_mode.tumor)
freq.normal <- gsub("CHR",1,freq_mode.normal)

t <- getdata(freq.tumor)
n <- getdata(freq.normal)

freq.common <- merge(t,n,by=c("Chr","Pos"),suffixes=c("T","N"))
ranges[with(ranges, startB <= startA & stopB >= stopA),]
res.test <- t(apply(freq.common[1:10,],1,function(x){myFETest(x[c(4,3)],x[c(7,6)])}))
rownames(res.test) <-freq.common$Pos[1:10]
chr1.test <- myChrFETest(13,1,snp.file)


res.test <- chr1.test

source("~/Dropbox/scripts/projNET/testCNV/CNVTest.r")
plotCNV22(FreqFileMode,snp.file,output)


###------------------------------testing

  grid.newpage()
  pushViewport(viewport(layout=grid.layout(6, 4)))
  par(mfrow=c(4,6),mai=c(0.1,0.2,0.1,0.1))
  for (i in 10:22) {
    # row <- floor((i - 1) / 4) + 1
    # column <- (1 + (i + 3) %% 4)
    freq.file <- gsub("CHR", i, FreqFileMode)
    print(getCNV(freq.file,snp.file)
      # , vp=vplayout(row, column)
      )
  }
  grid.text("test", vp=viewport(layout.pos.row = 6, layout.pos.col=3:4))
  dev.off()

}
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)














# tumor <- getSnpFreq22Data(freq_mode.tumor.file, snp.file)
# normal <- getSnpFreq22Data(freq_mode.normal.file, snp.file)