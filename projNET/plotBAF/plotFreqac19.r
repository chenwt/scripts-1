setwd("/Volumes/ys_lab_scratch/jh3283/net/AC19/res/")
require(ggplot2)
require(scales)

#####functions ####
getdata <- function(filename){
  Freqchr <- read.table(filename,skip=5,header=T,sep="\t")
  Freqchr <- Freqchr[Freqchr$Altreads
   > -1 & Freqchr$Altreads/Freqchr$Totalreads < 2,]
  Freqchr$Altfreq <- Freqchr$Altreads/Freqchr$Totalreads
  return(Freqchr)
}

freq_mode.file <- "../BAF/chrCHR.freq"
p_snp.file <- "~/net/p_snps.txt"

plotHistChr <- function(freq_mode.file,p_snp.file) {
###------------------------------ load individual snp
  p_snp <- read.table(p_snp.file,sep="\t",header=F)
  colnames(p_snp) <- c("Chr","Pos","Ale")
###------------------------------ load freq data
  
  Freq <- data.frame(Chr=integer(),Pos=integer(),Totalreads=integer(),Altreads=integer())
  for (i in 1:22){
    filename <-"../BAF/chrCHR.freq"
    filename <- gsub("CHR",i,filename)
    Freq <- rbind(Freq,getdata(filename))
  }
  Freqp_snp <- merge(x=Freq,y=p_snp,by=c("Chr","Pos"))
  Freqp_snp$Altfreq <- Freqp_snp$Altreads / Freqp_snp$Totalreads
  Freqp_snp$Zscore <- (Freqp_snp$Altfreq - mean(Freqp_snp$Altfreq))/sd(Freqp_snp$Altfreq)

  attach(Freqchr)
  qplot((Altfreq-mean(Altfreq))/sd(Altfreq),data=Freqchr,geom="histogram",binwidth=0.05)
  qplot(Altfreq,data=Freqchr,geom="histogram",binwidth=0.01)
  detach(Freqchr)
  
}


#### running in order
plotHist
pdf("ac19_Hist_22chr.pdf")
for (i in 1:22){
  data.plot <- Freqp_snp[Freqp_snp$Chr == i,]
#  attach(data.plot)
  plot1 <- qplot(Altfreq,data=data.plot,geom="histogram",binwidth=0.05)  + opts(title = paste("chr",i,sep=""))
  plot2 <- qplot(Zscore,data=data.plot,geom="histogram",binwidth=0.05)  + opts(title = paste("chr",i,sep=""))
  print(plot1)
  print(plot2)
  #detach(data.plot)
}
dev.off()

pdf("ac19_Point_22chr.pdf")
for (i in 1:22){
  data.plot <- Freqp_snp[Freqp_snp$Chr == i,]
  p <- ggplot(data=data.plot,aes(y=Altfreq,x=Pos)) + geom_point(col="blue") + opts(title=paste("chr",i,sep=""))
  p <- p + data.plo
  print(p)
}
dev.off()

reg.del <- read.table("../../AC3somatic_corrected_NoModLong-3-type.seg",sep="\t",header=T)
data.plot <- Freqp_snp[which(Freqp_snp$Chr == 3 & Freqp_snp$Pos >= 66643001 & Freqp_snp$Pos <=72869000),]
ggplot(data=data.plot,aes(y=Altfreq,x=Pos)) + geom_point(col="blue") + opts(title=paste("chr",3,sep=""))

for ( i in 1:nrow(reg.del)){
  data.plot <- Freqp_snp[which(Freqp_snp$Chr == reg.del$CHR[i] & Freqp_snp$Pos >= reg.del$BP1[i] & Freqp_snp$Pos <= reg.del$BP2[i]),]
  p <- ggplot(data=data.plot,aes(y=Altfreq,x=Pos)) + geom_point(col="blue") + opts(title=paste("chr",reg.del$CHR[i],sep=""))  
  p <- p + xlim(reg.del$BP1[i],reg.del$BP2[i])
  print(p)
}

####################
FreqtoBAFPlot <- function(freq.file, chr.num){
  freq <- read.table(freq.file, skip=5, header=TRUE)
  freq$Altfreq <- freq$Altreads / freq$Totalreads
  freq <- freq[freq$Altfreq > -1 & freq$Altfreq < 2, ]
  #plot <- ggplot(freq, aes(Pos, Altfreq)) + 
  #        geom_point(color=alpha("black", 1/250)) +
  #        labs(x=paste("Chr", chr.num), y=NULL) + 
  #        opts(plot.margin = unit(rep(0, 4), "lines")) 
  #plot <- ggplot(freq, aes(Pos, Altfreq)) + geom_point(color="#ffc0cb") +
  plot <- qplot(x=Pos,y=Altfreq,data=freq,geom="point", color="#ffc0cb", size=freq$Totalreads) + 
    scale_size(range = c(0.3, 3)) + 
    opts(panel.background = theme_rect(fill='white', color='black'))+
    #stat_binhex(bins=40) +
    #scale_fill_continuous(low="white", high="red") +
    labs(x=NULL, y=NULL) +
    opts(plot.margin = unit(rep(0, 4), "lines")) +
    opts(legend.position = "none") +
    opts(axis.text.x = theme_blank()) +
    opts(axis.text.y = theme_blank()) +
    opts(title=paste("Chr", chr.num)) +
    #scale_x_continuous(limits=c(0, 250000000), breaks=NA) +
    scale_y_continuous(breaks=NA)
  return(plot)
  #scale_x_continuous('x axis label')
}

FreqtoBAFPlot("../BAF/chr1.freq",1)
freq.file.model <- "../BAF/chrCHR.freq"
FreqBAFPlot22Chrs(freq.file.model,output.file="ac19_22chr_BAF.pdf",plot.title= "AC19 BAF")
