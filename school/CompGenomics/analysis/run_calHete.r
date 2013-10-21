#!/bin/Rscirpts
#Author: Jing He
#Date: Mar. 2013
#Last Update: Mar.24, 2013


# snp.file <-  "/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt"
# t.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC13/BAF/chr11.freq"
# n.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chr11.freq"
# bin.num <- 50
# install.packages("reshape2")
# require(reshape2)
# dcast(SNPdata, cut(AltFreq, breaks = c(0, 25, 100)) ~ status)

# x <- as.numeric(SNPdata$Altfreq )
# br <- seq(0, max(SNPdata$Altfreq), by=0.05)
# xt <- summary(cut(x[x!=0], br, labels=br[-1], include.lowest=T, ordered_result=T))

# vapply(SNPdata$Altfreq,mySNPHete,1)
# library(help="diptest")
# dip.test(SNPdata$Altfreq[SNPdata$Altfreq > 0])
# plot(bkde(v[61:70]),type="l",cex=3)



