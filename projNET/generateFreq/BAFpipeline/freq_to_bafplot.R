
FreqtoBAFPlot <- function(freq.file, chr.num){
  freq <- read.table(freq.file, skip=5, header=TRUE)
  freq$Altfreq <- freq$Altreads / freq$Totalreads
  freq <- freq[freq$Altfreq > 0 & freq$Altfreq < 1, ]
  #plot <- ggplot(freq, aes(Pos, Altfreq)) + 
  #        geom_point(color=alpha("black", 1/250)) +
  #        labs(x=paste("Chr", chr.num), y=NULL) + 
  #        opts(plot.margin = unit(rep(0, 4), "lines")) 
  plot <- ggplot(freq, aes(Pos, Altfreq)) + 
          opts(panel.background = theme_rect(fill='white', color='black'))+ 
          stat_binhex(bins=40) +
          scale_fill_continuous(low="white", high="red") +
          labs(x=NULL, y=NULL) +
          opts(plot.margin = unit(rep(0, 4), "lines")) +
          opts(legend.position = "none") + 
          opts(axis.text.x = theme_blank()) +
          opts(axis.text.y = theme_blank()) +
          opts(title=paste("Chr", chr.num)) +
          scale_x_continuous(limits=c(0, 250000000), breaks=NA) +
          scale_y_continuous(breaks=NA)
  return(plot)
  #scale_x_continuous('x axis label')
}

vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)

FreqBAFPlot22Chrs <- function(freq.file.model, output.file, plot.title){
#  freq.file.model <- "/ifs/scratch/c2b2/ngs_lab/db2175/BAF/DBSNPData/AC3DBSNPCommonSNPsfreqs/AC3dbSNP135commonchrCHR.freq"
  library(ggplot2)
  library(grid)
  library(scales)
  plot.list <- vector('list', 22)
  for(i in 1:22) {
    freq.file <- gsub("CHR", i, freq.file.model)
    plot.list[[i]] <- FreqtoBAFPlot(freq.file, i)
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



args <- commandArgs()
freq.file.model <- args[9]
output.file <- args[10]
plot.title <- args[11]
Routput.file <- args[12]
plot.title <- gsub("_", " ", plot.title)
sink(Routput.file)
FreqBAFPlot22Chrs(freq.file.model, output.file, plot.title)

