##!/usr/bin/Rscript
#Author: Jing He
#Date:Jul. 17th, 2013 
#Last Updated: Jul. 29rd, 2013
#Usage: 
#Description: Input: tumor acronym, prepare all compressed files and unzip them; output: a txt file with genes and Fisher's Method chi-square statistics
 

args <- commandArgs(TRUE)
if (is.null(args)){
  print("Please provide parameters")
  exit
}else{
  print(args)
}

sd <- "/ifs/scratch/c2b2/ac_lab/jh3283/projPanCancer/myethy/"
#tumAcro <- "brca"
tumAcro <- args[1]

cwd <- paste(sd,"/",tumAcro,sep="")
setwd(cwd)
fns <- dir(pattern=paste("tar.gz"))

# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")

cnvRead <- function(dpattern="analysis.hg19.seg", ofile, barcode=NULL) {
        require(GenomicFeatures)
        load("/ifs/scratch/c2b2/ac_lab/malvarez/databases/genomes/hg19-gene.rda")
        fn <- list.files(pattern=dpattern)
        if (!is.null(barcode)) fn <- sapply(colnames(dset), function(x, fn) fn[grep(x, fn)], fn=fn)
        if (length(fn)==0) stop("No files to process", call.=F)
        cnv <- NULL
        barcode <- NULL
  linkage <- list()
  pb <- txtProgressBar(max=length(fn), style=3)
        for (i in 1:length(fn)) {
                tmp <- strsplit(readLines(fn[i]), "\t")
                tmp <- lapply(tmp, function(x) {if(x[2]=="Y") {x<-NULL}; x}) # Totally cryptic anonymous function to please Mariano
                tmp <- tmp[-(which(sapply(tmp,is.null),arr.ind=TRUE))]
                if (length(tmp)==0) {next}
                cnvd <- lapply(1:length(tmp[[1]]), function(i, tmp) sapply(tmp, function(x, i) ifelse(i>2, as.numeric(x[i]), x[i]), i=i), tmp=
                names(cnvd) <- tmp[[1]]
                cnvd[[2]] <- paste("chr", cnvd[[2]], sep="")
                cnvd1 <- GRanges(seqnames=cnvd$chromosome, ranges=IRanges(start=cnvd$start, end=cnvd$stop), strand=rep("*", length(cnvd$chromo
                tmp <- as.matrix(findOverlaps(cnvd1, tx_by_gene))
    tmp[, 2] <- names(tx_by_gene)[tmp[, 2]]
    tmp <- cbind(tmp, cnvd$seg.mean[as.numeric(tmp[, 1])])
    tmp <- tmp[order(abs(as.numeric(tmp[, 3])), decreasing=T), ]
    tmp <- tmp[!duplicated(tmp[, 2]), ]
    linkage <- c(linkage, list(split(tmp[, 2], tmp[, 1])))
    genes <- unique(c(rownames(cnv), tmp[, 2]))
                cnv <- cbind(cnv[match(genes, rownames(cnv)), ], as.numeric(tmp[match(genes, tmp[, 2]), 3]))
                rownames(cnv) <- genes
                barcode <- c(barcode, cnvd[[1]][1])
                setTxtProgressBar(pb, i)
        }
        colnames(cnv) <- substr(barcode, 1, 16)
  names(linkage) <- colnames(cnv)
        save(cnv, linkage, file=paste(ofile, "-cnv.rda", sep=""))
}



getBetaParam <- function(x){ # x is a vector of values which to be fitted using beta distribution
    x <- as.numeric(as.character(na.omit(x)))
    mu <- mean(x)
    var <- var(x)
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))  
}
getNearestProb <- function(gene){
#res.temp <- data.frame("Composite.Element.REF"=factor(),"Beta_value"=integer(),"Gene_Symbol"=factor(),"Chromosome"=factor(),"Genomic_Coordinate"=factor(),"p_value"=integer())
  res.temp <- rep(NA,6)
  for (i in 1:length(gene)){
    g <- as.character(gene[i])
    genedata <- res[which(res$Gene_Symbol==g),]
    geneTss <- as.numeric(dataTSS$V2[which(dataTSS$V3==g)])
    if(dim(genedata)[1]>=1 & length(geneTss) >=1){
    res.temp <- rbind(res.temp,genedata[which.min(abs(sapply(genedata$Genomic_Coordinate,FUN=function(x){min(as.numeric(x) - as.numeric(geneTss))}))),])
   # res.temp[i,] <- genedata[which.min(abs(sapply(genedata$Genomic_Coordinate,FUN=function(x){min(as.numeric(x) - as.numeric(geneTss))}))),]
    }}
  colnames(res.temp) <- colnames(res)
  return(na.omit(res.temp))
}

getResult <- function(fn){
  subdir <- gsub(".tar.gz","",fn)
#  subdir <-  dir()[file.info(dir())$isdir]
  fn <- dir(path=subdir,pattern=paste("*",toupper(tumAcro),"*",sep=""))
  data <- na.omit(read.table(paste(subdir,"/",fn[1],sep=""),skip=1,sep="\t",header=T))
  data <- data[data$Gene_Symbol!="",]
  lambda <- mean(data$Beta)
  
  betaParam <- getBetaParam(data$Beta[data$Beta < 0.2])
  p_value <- pbeta(data$Beta_value, shape1=betaParam$alpha,shape2=betaParam$beta, lower.tail=F)
  res <- cbind(data,p_value=p_value)
  #z_score <- sapply(p_value, FUN=function(x){qnorm(x,lower.tail=FALSE)})
  res <- res[order(res$p_value),]
  res <- res[which(res$p_value <= 1e-4),]
  genes <- as.character(unlist(na.omit(unique(unlist(res$Gene_Symbol)))))
  #need further modification to integrate TSS information
  #dataTSS <- read.table("~/SCRATCH/ref/RefSeq_Gene_TSS.txt", sep="\t",skip=1,stringsAsFactors=F)
  genes <- as.character(unlist(na.omit(genes)))
  #result <- getNearestProb(genes)
  #write.table(result,paste("significant_probes_", tumAcro,".txt", sep=""),quote=F,sep="\t",row.names=F)
  #Just get the smallest p-value probe for one gene
  require(plyr)
  result <- arrange(res,p_value,desc(Beta_value))
  genes <- sort(unique(res$Gene_Symbol))
  result <- result[match(genes,result$Gene_Symbol),]
  print(paste("writing file...","significant_probes_topPvalue_", tumAcro,"_",substr(subdir,nchar(subdir)-8+1,nchar(subdir)),".txt",sep=""))
  write.table(result,paste("significant_probes_topPvalue_", tumAcro,"_",substr(subdir,nchar(subdir)-8+1,nchar(subdir)),".txt", sep=""),quote=F,sep="\t",row.names=F)
  #hist(qnorm(result$p_value,lower.tail=F),xlab="z-score",breaks=100,main="Z-socre Hist")
}

#write the files
sapply(fns,getResult)

#integrating meta studies
fns <- dir(pattern="significant_probes_topPvalue_brca")
dataMeta <- read.table(fns[1],header=T,sep="\t")

for (i in 2:length(fns)){
  if (i == 2){
    dataAll <- merge(dataMeta[,c(1,3,4,5,6)],read.table(fns[i],header=T,sep="\t")[,c(3,6)],by="Gene_Symbol") 
  }else{
    temp <- read.table(fns[i],header=T,sep="\t")[,c(3,6)]
    colnames(temp) <- c("Gene_Symbol",substr(fns[i],nchar(fns[i])-10+1,nchar(fns[i])-4))
    dataAll <- merge(dataAll,temp,by="Gene_Symbol")
  }
}

resultAll<- t(apply(dataAll,1,function(x){data <- sapply(x[5:length(x)],as.numeric);
					  Xsq <- round(-2 * sum(log(data)),digits=2);
					  res <- c(x[1:4],
					  fisher_Method=Xsq,
					  p_value= pchisq(Xsq,df=2*length(fns),lower.tail=F))
	      })) 
write.table(resultAll, paste(tumAcro,"_genelist_integrate_",ncol(dataAll) - 4, "_FisherM.txt",sep=""),sep="\t",quote=F,row.names=F)

