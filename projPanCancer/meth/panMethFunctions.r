##----------------------------
#functions for methylation


##----------------------------
#to Get best probe for each gene
#modified from cnvRead function in /ifs/scratch/c2b2/ac_lab/malvarez/panCancer/scripts/panCancerFunctions.r

methRead <- function(dpattern="HumanMethylation", ofile, barcode=NULL) {
#This function will save data to a rda file
# file header for 
# Composite.Element.REF Beta_value Gene_Symbol Chromosome Genomic_Coordinate
# cg00000029 0.05554673        RBL2         16           53468112
# cg00000108         NA     C3orf35          3           37459206
# cg00000109         NA      FNDC3B          3          171916037
        require(GenomicFeatures)
        load("/ifs/scratch/c2b2/ac_lab/malvarez/databases/genomes/hg19-gene.rda")
        fn <- list.files(pattern=dpattern)
        ##----------------------------
        fn <- fn[1:3]
        #choose those samples gived in barcode 
        if (!is.null(barcode)) fn <- sapply(colnames(dset), function(x, fn) fn[grep(x, fn)], fn=fn)
        if (length(fn)==0) stop("No files to process", call.=F)
        meth <- NULL
        barcode <- NULL
        linkage <- list()
        pb <- txtProgressBar(max=length(fn), style=3)
        for (i in 1:length(fn)) {
                tmp<- strsplit(readLines(fn[i]), "\t")
                barcode <- c(barcode, tmp[[1]][2])
                tmp <- tmp[2:length(tmp)]
                ## quality control
                tmp <- lapply(tmp, function(x) {if(x[4]=="Y") {x<-NULL}; x}) # get rid of chromosome Y meth
                tmp <- tmp[-(which(sapply(tmp,is.null),arr.ind=TRUE))] # get rid of records which has any null value
                if (length(tmp)==0) {next}
                methd <- lapply(1:length(tmp[[1]]), function(i, tmp) sapply(tmp, function(x, i) ifelse(i!=1 & i!=3, as.numeric(x[i]), x[i]), i=i), tmp=tmp[-1])
                names(methd) <- tmp[[1]]
                methd[[4]] <- paste("chr", methd[[4]], sep="")
                methd1 <- GRanges(seqnames=methd$Chromosome, ranges=IRanges(start=methd$Genomic_Coordinate - 1000, end=methd$Genomic_Coordinate + 1), strand=rep("*", length(methd$Chromosome)))
                tmp <- as.matrix(findOverlaps(methd1, tx_by_gene))
                tmp[, 2] <- names(tx_by_gene)[tmp[, 2]]
                tmp <- cbind(tmp, methd$Beta_value[as.numeric(tmp[, 1])])
                tmp <- tmp[order(as.numeric(tmp[, 3]), decreasing=T), ]
                tmp <- tmp[!duplicated(tmp[, 2]), ]
                linkage <- c(linkage, list(split(tmp[, 2], tmp[, 1])))
                genes <- unique(c(rownames(meth), tmp[, 2]))
                meth <- cbind(meth[match(genes, rownames(meth)), ], as.numeric(tmp[match(genes, tmp[, 2]), 3]))
                rownames(meth) <- genes
                setTxtProgressBar(pb, i)
        }
        colnames(meth) <- substr(barcode, 1, 16)
        names(linkage) <- colnames(meth)
        save(meth, linkage, file=paste(ofile, "-meth.rda", sep=""))
}