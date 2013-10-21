#!/usr/bin/Rscript
#Author: Jing He
#Date: Jun 15,2013 
#Last Updated: Jun 16,2013
#Usage: Rscript getMAF_aml.R <pid-sampleid> <NumCHR> 
#Description: generate frequency file using pileup files
# this code is based on codes from Daniel's BAFCodeExtrats.r in ~/scripts/projNET/BAFCodeExtracts.J.r

 ###--------------------------------------
 ## get arguments from command line
args <- commandArgs(TRUE)
if (is.null(args)){
  print("Please provide parameters")
  exit
}else{
  print(args)
}
 ###--------------------------------------

##--------------------------------------
#get variables 
pid <- args[1]
chr <- args[2]

##--------------------------------------
#set up environment variable
Sys.setenv(REF="/ifs/scratch/c2b2/ac_lab/jh3283/ref/GRCh37-lite.fa")
Sys.setenv(DBSNPMODEL="/ifs/scratch/c2b2/ac_lab/jh3283/ref/DBSNPCommonSNPsbedfiles/dbSNP135commonchrCHR.bed")

##--------------------------------------
# main procedure
get_freq <- function(pid,chr){

  pileup.df <- read.delim(paste(pid,"_",chr,".pileup",sep=""),sep="\t",header=FALSE)
  colnames(pileup.df) <- c("RNAME", "POS", "RBASE", "NUMREADS", "BASEREADS", "BASEQUALITIES")
 for (i in 1:dim(pileup.df)[1]){
      pileup.df$NUMMINORALLELE[i] <- sort(unlist(
        PileupToCounts(pileup.df$BASEREADS[i], pileup.df$RBASE[i], 
                       pileup.df$NUMREADS[i])), decreasing=TRUE)[2]
    }
    freq.df <- pileup.df[,c(1,2,3,4,7)] 
    write.table(freq.df,paste("PAEEYP-04A_chr",chr,".freq",sep=""),row.names=F,quote=F,sep="\t")
}


##--------------------------------------
#functions 

##--------------------------------------
# prepare DBSNP BEDfile for pileup
##--------------------------------------
#samtools pileup

##--------------------------------------
# pileup to counts
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

get_freq(pid,chr)
