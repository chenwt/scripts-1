Sys.setenv(REF="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa")
Sys.setenv(DBSNPMODEL="/ifs/scratch/c2b2/ac_lab/jh3283/ref/DBSNPCommonSNPsbedfiles/dbSNP135commonCHR.txt")
bam.file <- "/ifs/home/c2b2/ac_lab/jh3283/SCRATCH/projAML/WXS/callVars/PAEEYP-04A.rmdup.new.bam"

test <- function(bam.file){

  DBSNPs.df <- GetDBSNPs(1, 100000, 200000)
  pileup.df <- RangesDFToPileupDF(bam.file, DBSNPs.df)    

}


PlotSNPs <- function(chromosome, start, end){
	  # get all DBSNPs between chromosome, start and end
	  DBSNPs.df <- GetDBSNPs(chromosome, start, end) 
	  # get a pileup for each DBSNP
	  pileup.df <- RangesDFToPileupDF(bam.file, DBSNPs.df)    
	  # for debugging purposes, save the pileup so don't have to generate it anew each time
	  save(pileup.df, file="pileup.temp")
	  #load("pileup.temp")
	  pileup.df$BASEREADS <- as.character(pileup.df$BASEREADS)
	  #TO DO (low priority):  make sure not counting N as minor allele
	  # calculate minor allele frequency for each position
	  for (i in 1:dim(pileup.df)[1]){
	    #print(pileup.df$BASEREADS[i])
	    pileup.df$NUMMINORALLELE[i] <- sort(unlist(
	      PileupToCounts(pileup.df$BASEREADS[i], pileup.df$RBASE[i], 
	                     pileup.df$NUMREADS[i])), decreasing=TRUE)[2]
	  }
	  # plot the SNPs
	  par(fig=c(0, 1, 0.7, 1), new=TRUE)
	  par(mar=c(0, 0.5, 2, 0.5))
	  SNPs.to.plot <- pileup.df[pileup.df$NUMREADS > 10 & 
	                              pileup.df$NUMMINORALLELE < pileup.df$NUMREADS &
	                              pileup.df$NUMMINORALLELE > 0, ]
	  plot(SNPs.to.plot$POS, SNPs.to.plot$NUMMINORALLELE / 
	         SNPs.to.plot$NUMREADS, yaxt='n', xaxt='n',
	       xlim=c(start, end), ylim=c(0, 1), col="red")
	  title(main=paste("chr", chromosome, ":", start, "to", end))
	#browser()
}


GetDBSNPs <- function(chromosome,start, end){
	
	DBSNP.file.model <- Sys.getenv("DBSNPMODEL")
  DBSNP.file <- gsub("CHR", chromosome, DBSNP.file.model)
  	col.names <- c("CHROM", "CHROMSTART", "CHROMEND", "REFUCSC", "AVHET")
  	
  	command <- paste("cat ", DBSNP.file, " |awk '$2 >= ", format(start, scientific=FALSE),
  					" && $2 <= ", format(end, scientific=FALSE), " { print $0}'", " | cut -f1,2,3,4,5", sep="")
  	DBSNP.pipe <- pipe(command)
  	dbsnp.tmp <- readLines(DBSNP.pipe)
  	close(DBSNP.pipe)

  	if(!length(dbsnp.tmp)){
  		reads.df <- data.frame(CHROM=character(), CHROMSTART=integer(), CHROMEND=integer(), 
                           		REFUCSC=character(),AVHET=character())
  	}else {
  		    reads.df <- read.delim(file=textConnection(dbsnp.tmp), sep="\t",col.names=col.names, header=FALSE) 
  	}
  	return(reads.df)
}

RangesDFToPileupDF <- function(bam.file, ranges.df){
	ranges.df$CHROM <- gsub("chr", "", ranges.df$CHROM)
	write.table(ranges.df[, c("CHROM", "CHROMSTART", "CHROMEND")], "temp.bed", quote=FALSE, sep="\t", row.names=FALSE, 
	              col.names=FALSE)
	max.depth <- 10000
	min.map.q <- 0
	min.base.q <- 0
	#use environment variable
	ref <- Sys.getenv("REF")
	command <- paste("samtools mpileup -DS -C 5 -d ", max.depth, " -q ", min.map.q, " -Q ", min.base.q, " -f ", ref, 
	                   " -l temp.bed ", bam.file, sep="")

	#print(command)
	#add support for empty range
	#to work in read.table, must set quote="", comment.char="" and potentially 
	#strip.white=TRUE
	pileup.df <- read.delim(pipe(command), sep="\t", header=FALSE)
	colnames(pileup.df) <- c("RNAME", "POS", "RBASE", "NUMREADS", "BASEREADS", "BASEQUALITIES")
	
	return(pileup.df)
}



PileupToCounts <- function(base.reads, ref.base, num.reads){
  #find indels
  # browser()
  indel.matches <- gregexpr("(-|\\+)\\d+", base.reads)[[1]]
  if (!indel.matches[1] == -1){    # are there any indels?
    indel.match.lengths <- attr(indel.matches, "match.length")
    num.matches <- length(indel.matches)
    indel.lengths <- rep(0, num.matches)
    for (i in 1:num.matches){
      indel.lengths[i] <- as.integer(substr(base.reads, 
                                            indel.matches[i] + 1, 
                                            indel.matches[i] + 1 + 
                                              indel.match.lengths[i] - 2))
    }
    #preserve reads before the first indel
    base.reads.without.indels <- substr(base.reads, 1, indel.matches[1] - 1)
    if (num.matches > 1){
      for (i in 2:num.matches){
        base.reads.without.indels <- paste(base.reads.without.indels, 
                                           substr(base.reads, indel.matches[i - 1] + 
                                                    indel.match.lengths[i - 1] + 
                                                    indel.lengths[i - 1], 
                                                  indel.matches[i] - 1), sep="")
      }
    }
    #preserve reads after the last indel
    base.reads.without.indels <- paste(base.reads.without.indels, 
                                       substr(base.reads, indel.matches[num.matches] +
                                                indel.match.lengths[num.matches] +
                                                indel.lengths[num.matches], 
                                              nchar(as.character(base.reads))), sep="")
  } else base.reads.without.indels <- base.reads
  #remove confounding characters associated with mapping qualities
  base.reads.without.indels <- 
    gsub("\\^\\S", "", base.reads.without.indels)      
  counts.list <- vector('list', 5)
  base.names <- c("A", "C", "G", "N", "T")
  names(counts.list) <- base.names
  total.count <- 0
  for (j in base.names){
    reg.ex <- paste("[^", j, tolower(j), "]", sep="")
    count <- nchar(gsub(reg.ex, "", base.reads.without.indels))
    counts.list[[j]] <- count
    total.count <- total.count + count 
  }
  counts.list[[toupper(ref.base)]] <- num.reads - total.count 
  return(counts.list)
}