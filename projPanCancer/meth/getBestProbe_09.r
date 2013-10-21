##!/usr/bin/Rscript
#Author: Jing He
#Date:Jul. 17th, 2013 
#Last Updated: Jul. 29rd, 2013
#Usage: 
#Description: Input: tumor acronym, prepare all compressed files and unzip them; output: a txt file with genes and Fisher's Method chi-square statistics
 
sd <- "/ifs/scratch/c2b2/ac_lab/jh3283/projPanCancer/meth/"
tumAcro <- "coad"

cwd <- paste(sd,"/",tumAcro,sep="")
setwd(cwd)

# system("for f in $(ls jhu-usc.edu_COAD.HumanMethylation*tar.gz) ; do tar -xvf ${f} -C . ; done")

# subfns <- system("tar -tvf jhu-usc.edu_COAD.HumanMethylation27.Level_3.1.3.0.tar.gz | cut -d' ' -f7| cut -d/ -f2|grep HumanMethylation",
		# intern = TRUE)
subd <- gsub(".tar.gz","",list.files(pattern="tar.gz"))


# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
fnexp <- paste("~/SCRATCH/projPanCancer/PANCANCER/",tumAcro,"/",tumAcro,"_tcga_rnaseq.exp",sep="")
dset <- read.table(fnexp, header=T)


methRead()
ofile <- tumAcro
for (i in 1:length(subd))
{
	subdf <- paste(cwd,"/",subd[i],sep="")
	setwd(file.path(cwd, subd[i]))
	subdfs <- dir(path=subdf, pattern=toupper(tumAcro))
	dpattern <- "HumanMethylation"
	
	barcode=NULL
	if(i == 1) {
		meth <- methRead(dpattern=dpattern,barcode=NULL)$meth
		# linkage <- methRead(dpattern=dpattern,barcode=NULL)$linkage

	}
	
	tmp <- methRead(dpattern=dpattern,barcode=NULL)$meth
	genes <- unique(c(rownames(meth), rownames(tmp)))
	meth <- cbind(meth[match(genes, rownames(meth)),],tmp[match(genes, rownames(meth)),])
	save(meth, linkage, file=paste(ofile, "-meth.rda", sep=""))
	# linkage <- methRead(dpattern=dpattern,barcode=NULL)$linkage
}

methRead <- function(dpattern=dpattern, barcode=NULL) {
        require(GenomicFeatures)
        load("/ifs/scratch/c2b2/ac_lab/malvarez/databases/genomes/hg19-gene.rda")
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
                tmp <- lapply(tmp, function(x) {if(x[4]=="Y") {x<-NULL}; x}) # Totally cryptic anonymous function to please Mariano
                tmp <- tmp[-(which(sapply(tmp,is.null),arr.ind=TRUE))]
                if (length(tmp)==0) {next}
                methd <- lapply(1:length(tmp[[1]]), 
                	function(i, tmp) sapply(tmp[2:length(tmp)], 
                		function(x, i) ifelse(i>1 & i!=3, as.numeric(x[i]), x[i]), i=i), tmp=tmp[-1])
                names(methd) <- tolower(tmp[[2]])
                methd[[4]] <- paste("chr", methd[[4]], sep="")
                methd1 <- GRanges(seqnames=methd$chromosome, 
                			ranges=IRanges(start=(methd$genomic_coordinate - 1000), 
                							end=(methd$genomic_coordinate + 1)),  
                			strand=rep("*", length(methd$chromosome)))
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




