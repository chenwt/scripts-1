require(plyr)
require(ggplot2)
setwd("/Volumes/ys_lab_scratch/jh3283/school/compGenomic/cgquery")
source("~/Dropbox/scripts/myUseful/R/myIO.r")
fns <- c(dir()[grep("*exome.txt",dir())],dir()[grep("*WXS.txt",dir())])
fns <- c(fns,dir()[grep("*RNAseq.txt",dir())])
data_df <- ldply(fns, .fun=function(x){read.table(x,header=T,stringsAsFactors=F,sep="\t")})
plot_cnt <- table(data_df[,c("library_strategy","disease_abbr")])
# plot(plot_cnt,type="h",col=c(5:7),main="Sample Size")

# ggplot(as.data.frame.matrix(plot_cnt),aes(library_strategy,Freq,color=disease_abbr)) + geom_hist()

mosaicplot(~  disease_abbr + library_strategy, data = plot_cnt,col=c(5,7),main="Sample Size") 
text( .4,.85, "339",col = "red", cex = 1) 
text( .7,.85, "173 ",col = "red", cex = 1) 
text( .92,.85, "54 ",col = "red", cex = 1) 
text( .4,.3, "1000",col = "blue", cex = 1) 
text( .7,.3, "365",col = "blue", cex = 1) 
text( .92,.3, "219",col = "blue",  cex = 1) 

###------------------------------ plot tissue type
code <- read.table("../codeTable/sampleType.txt",sep="\t",header = T)
# head(code)


data_df <- merge(x=data_df,y=code, by.x="sample_type",by.y="Code")
plot_cnt <- table(data_df[,c("disease_abbr","Short.Letter.Code")])
plot_cnt <- as.data.frame.matrix(plot_cnt)

library("ggplot2")
require(grid)
ylab <- lapply(1:3,FUN=function(x){names(plot_cnt[,which(plot_cnt[x,] != 0)])})
par(mfrow=c(3,1))
plot(plot_cnt,col=c(4:8),main="Sample type")



	