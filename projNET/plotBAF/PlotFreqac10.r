setwd("/Volumes/ys_lab_scratch/jh3283/net/AC10/res/")
getdata <- function(filename){
  Freqchr <- read.table(filename,skip=5,header=T,sep="\t")
  Freqchr <- Freqchr[Freqchr$Altreads > -1 & Freqchr$Altreads/Freqchr$Totalreads < 2,]
  Freqchr$Altfreq <- Freqchr$Altreads/Freqchr$Totalreads
  return(Freqchr)
}

Freq <- data.frame(Chr=integer(),Pos=integer(),Totalreads=integer(),Altreads=integer())
for (i in 1:22){
  filename <-"../BAF/chrCHR.freq"
  filename <- gsub("CHR",i,filename)
  Freq <- rbind(Freq,getdata(filename))
}


ac1snp <- data.frame(Chr=integer(),Pos=integer(),Alt=character())
ac1snp <- read.table("~/net/AC1SNPs.txt",sep="\t",header=F)
colnames(ac1snp) <- c("Chr","Pos","Ale")

FreqAc1Snp <- merge(x=Freq,y=ac1snp,by=c("Chr","Pos"))
FreqAc1Snp$Altfreq <- FreqAc1Snp$Altreads / FreqAc1Snp$Totalreads
FreqAc1Snp$Zscore <- (FreqAc1Snp$Altfreq - mean(FreqAc1Snp$Altfreq))/sd(FreqAc1Snp$Altfreq)


pdf("ac10_Hist_22chr.pdf")
for (i in 1:22){
  data.plot <- FreqAc1Snp[FreqAc1Snp$Chr == i,]
  #  attach(data.plot)
  plot1 <- qplot(Altfreq,data=data.plot,geom="histogram",binwidth=0.05)  + opts(title = paste("chr",i,sep=""))
  plot2 <- qplot(Zscore,data=data.plot,geom="histogram",binwidth=0.05)  + opts(title = paste("chr",i,sep=""))
  print(plot1)
  print(plot2)
  #detach(data.plot)
}
dev.off()

pdf("ac10_Point_22chr.pdf")
for (i in 1:22){
  data.plot <- FreqAc1Snp[FreqAc1Snp$Chr == i,]
  p <- ggplot(data=data.plot,aes(y=Altfreq,x=Pos)) + geom_point(col="blue") + opts(title=paste("chr",i,sep=""))
  print(p)
}
dev.off()


reg.del <- read.table("../../AC3somatic_corrected_NoModLong-3-type.seg",sep="\t",header=T)

pdf("ac10_delreg_snp.pdf")
for ( i in 1:nrow(reg.del)){
  data.plot <- FreqAc1Snp[which(FreqAc1Snp$Chr == reg.del$CHR[i] & FreqAc1Snp$Pos >= reg.del$BP1[i] & FreqAc1Snp$Pos <= reg.del$BP2[i]),]
  p <- ggplot(data=data.plot,aes(y=Altfreq,x=Pos)) + geom_point(col="blue") + opts(title=paste("chr",reg.del$CHR[i],sep=""))  
  p <- p + xlim(reg.del$BP1[i],reg.del$BP2[i])
  print(p)
}
dev.off()


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
    scale_size(range = c(0.2, 3)) + 
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

FreqtoBAFPlot("../BAF/chr10.freq",10)
freq.file.model <- "../BAF/chrCHR.freq"
FreqBAFPlot22Chrs(freq.file.model,output.file="ac10_22chr_BAF.pdf",plot.title= "AC10 AltFrequent")

