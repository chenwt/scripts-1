##!/usr/bin/Rscript
#Author: Jing He
#Date:Jul. 17th, 2013 
#Last Updated: Sep 05, 2013
#Usage: 
#Description: Input: tumor acronym, prepare all compressed files and unzip them; 
#               output: a txt file with genes and Fisher's Method chi-square statistics
 
# TODO:IMPROVE the effeicence!!! make this fast!
##----------------------------
# library(BSgenome)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(GenomicFeatures)

methRead <- function(dpattern=dpattern, barcode=NULL) {
        # require(GenomicFeatures)
        geneAnno <- read.table("~/SCRATCH/ref/aclab/genes_codon.tab",header=T,sep="\t",stringsAsFactors=F)
        load("~/SCRATCH/ref/aclab/hg19-gene.rda")
        fn <- list.files(pattern=dpattern)
        # Hybridization REF       TCGA-A6-2672-11A-01D-1551-05    TCGA-A6-2672-11A-01D-1551-05    TCGA-A6-2672-11A-01D-1551-05    TCGA-A6-2672-11A-01D-1551-05
    		# Composite Element REF   Beta_value      Gene_Symbol     Chromosome      Genomic_Coordinate
    		# cg00000292      0.511852232819811       ATP2A1  16      28890100
    		# barcode <- gsub("-",".",strtrim(x=unlist(strsplit(fn[1],split="\\."))[6],width=16))
    
        if (!is.null(barcode)) fn <- sapply(colnames(dset), function(x, fn) fn[grep(x, fn)], fn=fn)
        if (length(fn)==0) stop("No files to process", call.=F)
        meth <- NULL
        barcode <- NULL
  		  linkage <- list()
  	   	pb <- txtProgressBar(max=length(fn), style=3)
        for (i in 1:length(fn)) {
                tmp <- strsplit(readLines(fn[i]), "\t")
                tmp <- lapply(tmp, function(x) {if(x[4]=="Y") {x<-NULL}; x}) # exclude Y chromosome
                tmp <- tmp[-(which(sapply(tmp,is.null),arr.ind=TRUE))]

                if (length(tmp)==0) {next}
                methd <- lapply(1:length(tmp[[1]]), 
                            	function(i, tmp) sapply(tmp[2:length(tmp)], 
                              		function(x, i) ifelse(i>1 & i!=3, as.numeric(x[i]), x[i]), i=i), tmp=tmp[-1])
                names(methd) <- tolower(tmp[[2]])
                methd[[4]] <- paste("chr", methd[[4]], sep="")
                methd <- lapply(methd,function(x) x[-which(methd$gene_symbol %in% geneAnno$gene_name ==FALSE,arr.ind=T)]) 
                methd$strand <- geneAnno$strand[which(methd$gene_symbol %in% geneAnno$gene_name == TRUE,arr.ind=T)]
                
                methd$start <-   methd$genomic_coordinate - 2
               

                methd$d2TSS[which(methd$strand == "+",arr.ind=T)] <- 
                    methd$genomic_coordinate[which(methd$strand == "+",arr.ind=T)] - geneAnno$start[which(methd$strand == "+",arr.ind=T)]
                methd$d2TSS[which(methd$strand == "-",arr.ind=T)] <- 
                    geneAnno$start[which(methd$strand == "-",arr.ind=T)] - methd$genomic_coordinate[which(methd$strand == "-",arr.ind=T)]

                  # TODO make the methd2 dataframe into a data frame
                methd2 <- data.frame(t(sapply(unique(methd$gene_symbol),function(x){
                                   tmp_idx <- which(methd$gene_symbol==x) 

                                   tmp <- data.frame(gene_symbol=methd$gene_symbol[tmp_idx],chromosome=methd$chromosome[tmp_idx],
                                                 start=as.integer(methd$start[tmp_idx]),genomic_coordinate=as.integer(methd$genomic_coordinate[tmp_idx]) + 1 ,
                                                 beta_value=as.numeric(methd$beta_value[tmp_idx]),d2TSS=as.integer(methd$d2TSS[tmp_idx]),
                                                 strand=methd$strand[tmp_idx])

                                    if(length(which(tmp[,6] <0 &tmp[,6] > -1000,arr.ind=T))){
                                      result <- tmp[which.min(tmp[which(tmp[,6] < 0 & tmp[,6] > -1000,arr.ind=T),]),]
                                    }else{
                                      result <- tmp[which.min(abs(tmp[,6])),]
                                    }
                                    return(result)
                              })))
               
                methd_temp <- list(gene_symbol=unlist(methd2[,1]),chromosome=unlist(methd2[,2]),
                              start=unlist(methd2[,3]),end=unlist(methd2[,4]),
                              beta_value=unlist(methd2[,5]), d2TSS=unlist(methd2[,6]),
                              strand=unlist(methd2[,7]))
                methd2 <- methd_temp

                methd1 <- GRanges(seqnames=methd2$chromosome, 
                            			ranges=IRanges(start=methd2$start, 
                                    							end=methd2$end),  
                            			strand=methd2$strand)

                tmp <- as.matrix(findOverlaps(methd1, tx_by_gene))
      			    tmp[, 2] <- names(tx_by_gene)[tmp[,2]]
      			    tmp <- cbind(tmp, methd$beta_value[as.numeric(tmp[, 1])])
      			    tmp <- tmp[order(abs(as.numeric(tmp[, 3])), decreasing=T), ]
      			    tmp <- tmp[!duplicated(tmp[, 2]), ]
      			    linkage <- c(linkage, list(split(tmp[, 2], tmp[, 1])))
      			    genes <- unique(c(rownames(meth), tmp[, 2]))
                meth <- cbind(meth[match(genes, rownames(meth)), ], as.numeric(tmp[match(genes, tmp[, 2]), 3]))
                rownames(meth) <- genes
                barcode <- c(barcode, gsub("-",".",unlist(strsplit(fn[i],"\\."))[6]))
                setTxtProgressBar(pb, i)
        }
        colnames(meth) <- substr(barcode, 1, 16)
  		names(linkage) <- colnames(meth)
  		return(list(meth=meth,linkage=linkage))
        # save(meth, linkage, file=paste(ofile, "-meth.rda", sep=""))
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

