# setwd("/Volumes/ys_lab_scratch/jh3283/net/AC19/res/")
require(ggplot2)
require(scales)

# freq_mode.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC19/BAF/chrCHR.freq"
# p_snp.file <- "/Volumes/ys_lab_scratch/jh3283/net/AC1SNPs.txt"
#####functions ####
getdata <- function(filename){
  Freqchr <- read.table(filename,skip=5,header=T,sep="\t")
  Freqchr <- Freqchr[Freqchr$Altreads > -1 & Freqchr$Altreads/Freqchr$Totalreads < 2,]
  Freqchr$Altfreq <- Freqchr$Altreads/Freqchr$Totalreads
  return(Freqchr)
}


plotHist22Chr <- function(freq_mode.file, p_snp.file, out.name) {
###------------------------------ load individual snp
  p_snp <- read.table(p_snp.file,sep="\t",header=F)
  colnames(p_snp) <- c("Chr","Pos","Ale")
###------------------------------ load 22 freq data
  
  Freq <- data.frame(Chr=integer(),Pos=integer(),Totalreads=integer(),Altreads=integer())
  for (i in 1:22){
    filename <- freq_mode.file
    filename <- gsub("CHR",i,filename)
    Freq <- rbind(Freq,getdata(filename))
  }
  Freqchr <- merge(x=Freq,y=p_snp,by=c("Chr","Pos"))
  Freqchr$Altfreq <- Freqchr$Altreads / Freqchr$Totalreads
  Freqchr$Zscore <- (Freqchr$Altfreq - mean(Freqchr$Altfreq))/sd(Freqchr$Altfreq)

  pdf(paste("Hist_",out.name,"_22Chr.pdf",sep=""))
  for (i in 1:22){
    data.plot <- Freqchr[Freqchr$Chr == i,]
    # attach(data.plot)
    plot1 <- qplot(Altfreq,data=data.plot,geom="histogram",binwidth=0.05)  + opts(title = paste("chr",i,sep=""))
    plot2 <- qplot(Zscore,data=data.plot,geom="histogram",binwidth=0.05)  + opts(title = paste("chr",i,sep=""))
    print(plot1)
    print(plot2)
    #detach(data.plot)
  }
  dev.off()
}


FreqtoBAFPlot <- function(freq.file,snp.file ,chr.num){
  freq <- read.table(freq.file, skip=5, header=TRUE)
  freq$Altfreq <- freq$Altreads / freq$Totalreads
  freq <- freq[freq$Altfreq > -1 & freq$Altfreq < 2, ]

  print(snp.file)
  snp <- read.table(snp.file,sep="\t",header=F)
  colnames(snp) <- c("Chr","Pos","Ale")
###------------------------------ load 22 freq data
  
  snpFreq <- merge(x=freq,y=snp,by=c("Chr","Pos"))
  snpFreq$Altfreq <- snpFreq$Altreads / snpFreq$Totalreads
  # snpFreq$Zscore <- (Freqchr$Altfreq - mean(Freqchr$Altfreq))/sd(Freqchr$Altfreq)

  plot <- qplot(x=Pos,y=Altfreq,data=snpFreq,geom="point", color="#ffc0cb"
    , size=snpFreq$Totalreads
    # , size=1
    ) + 
    scale_size(range = c(0.5, 3)) + 
    opts(panel.background = theme_rect(fill='white', color='black'))+
    labs(x=NULL, y=NULL) +
    opts(plot.margin = unit(rep(0, 4), "lines")) +
    opts(legend.position = "none") +
    opts(axis.text.x = theme_blank()) +
    opts(axis.text.y = theme_blank()) +
    opts(title=paste("Chr", chr.num)) +
    scale_y_continuous(breaks=NA)
  return(plot)
}


FreqBAFPlot22Chrs <- function(freq.file.model, snp.file, output.file, plot.title){
  library(ggplot2)
  library(grid)
  library(scales)
  plot.list <- vector('list', 22)
  for(i in 1:22) {
    freq.file <- gsub("CHR", i, freq.file.model)
    plot.list[[i]] <- FreqtoBAFPlot(freq.file,snp.file, i)
  }
  pdf(output.file)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(6, 4)))
  for (i in 1:22) {
    row <- floor((i - 1) / 4) + 1
    column <- (1 + (i + 3) %% 4)
    print(plot.list[[i]], vp=vplayout(row, column))
  }
  grid.text(plot.title, vp=viewport(layout.pos.row = 6, layout.pos.col=3:4))
  dev.off()
}


vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)

