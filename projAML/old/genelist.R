setwd("c:/Documents and Settings/jh3283/My Documents/Dropbox/1Research/AML/")
# setwd("c:/Users/Sam Ni/Dropbox/1Research/AML/data/")

gl<- read.table("/BM.AML.gl.txt",skip=1,sep = "\t")
head(gl)
colnames(gl) <- c("probe.id","zf","pf","zd","pd","gene")
gltop200 <- gl[order(gl$pd)[1:200],c(5,6)]

require(ggplot2)
ggplot(gl,aes(x=pd)) + geom_histogram()

	Get_topGL <- function(gl,cutoff = 200)
	{
		gltop <- list()
		colnames(gl) <- c("probe.id","zf","pf","zd","pd","gene")
		gltop<- gl[order(gl$pd)[1:cutoff],6]
		return(gltop)
	}

tmp <- dir()
fns <- tmp[grep("gl",list.files())]
res.gl <- data.frame(rep(NA,200))


for (i in 1:length(fns)) {
	gl<- read.table(fns[i],skip=1,sep = "\t")
	mygl <- Get_topGL(gl,200)	
	res.gl[,i] <- mygl
}

colnames(res.gl) <- gsub("txt","",fns)
res.gl

write.table(res.gl,"Stirward_5_HUGO.gmx",row.names=F,sep= "\t")

