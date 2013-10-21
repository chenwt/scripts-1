##!/usr/bin/Rscript
#Author: Jing He
#Date: 
#Last Updated:
#Usage:
#Description:

# FindGenesInCNVs.4.19.2013.R takes a segmentation file
# and writes a list of the genes that overlap with any of the CNVs
# in the segmentation file

# it takes three arguments:

# 1. segmentation file (four tab-delimited columns CHR\tBP1\tBP2\tTYPE) 
#    where TYPE is the copy number estimate
# 2. ref seq database downloaded from UCSC 
# (3rd column CHR, 5th column genestartpos, 6th column geneendpos, 
# 13th column genename)
#     see /ifs/scratch/c2b2/ngs_lab/db2175/Projects/PersonalizedCancer/Ref/refGene.3.10.2013.txt
#     for an example
# 3. output file to write result to (writes 5 tab-delimited columns, 
#    one for each gene in a CNV, 
#    CHR\tCNVTYPE\tCNVSTARTPOS\tCNVENDPOS\tGENENAME)

# The script will output a summary like this:
#  2013-04-19 17:00:45 Script FindGenesInCNVs.4.19.2013.R: Marking genes 
#  overlapping with the cnvs listed in/ifs/scratch/c2b2/ngs_lab/db2175
#  /Projects/PersonalizedCancer/Data/regions.final.txt using the genes 
#  listed in /ifs/scratch/c2b2/ngs_lab/db2175/Projects
#  /PersonalizedCancer/Ref/refGene.3.10.2013.txt (table 
#  refGene (hg19/Genes and Gene Prediction Tracks/RefSeq Genes/refGene) 
#  downloaded on 3/10/2013 from genome.ucsc.edu).  Writing output to 
#  /ifs/scratch/c2b2/ngs_lab/db2175/Projects/PersonalizedCancer/Data
#  /regions.final.genes.txt


GetChromosomeNumber <- function(chromosome){
    # chromosome fields in the refseq.file can be, e.g., chr6_apd_hap1
    # (mostly on chromosome 6)
    gsub("chr", "", unlist(strsplit(chromosome, "_"))[1])
}

FindGenesInCNVs <- function(intervals.file, refseq.file,
			                                genes.in.cnvs.file){
    #refseq.readme.first.line <- readLines(
    #  paste(refseq.file, ".README", sep=""), 
    #  n=1)
    library(lubridate)
  #cat(as.character(now()), " Script MarkGenesInCNVs.4.19.2013.R: Marking genes ", 
      #"overlapping with the cnvs listed in ", 
      #intervals.file, " using the genes listed in ", refseq.file, " (", 
     # refseq.readme.first.line, ").  Writing output to ", 
      #genes.in.cnvs.file, "\n", sep="")
  intervals.df <- read.table(intervals.file, header=T)
    # column 13 of the refseq.file is name 2 (Alternate name, e.g., 
    # gene_id from GTF)
    refseq.df <- read.table(refseq.file, 
	          stringsAsFactors=F)[, c(1, 2, 3, 4)]
    colnames(refseq.df) <- c("CHR", "GENESTARTPOS", "GENEENDPOS", 
			                                "GENENAME")
      #refseq.df$CHR <- sapply(refseq.df$CHR, GetChromosomeNumber)
      #intervals.df$INTERVAL <- paste(intervals.df$CHR, ":", 
      #                               intervals.df$BP1, "-",
      #                               intervals.df$BP2, sep="")
      num.intervals <- dim(intervals.df)[1]
      cat("CHR\tCNVTYPE\tCNVSTARTPOS\tCNVENDPOS\tGENENAME\n", 
	        file=genes.in.cnvs.file)
        for (i in 1:num.intervals){
	      chr <- intervals.df$CHR[i]
	    cnv.start.pos <- intervals.df$BP1[i]
	        cnv.end.pos <- intervals.df$BP2[i]
	        genes.in.interval.df <- refseq.df[refseq.df$CHR == chr &
				                  refseq.df$GENESTARTPOS <= cnv.end.pos & 
						 refseq.df$GENEENDPOS >= cnv.start.pos, ]
		    if (length(unique(genes.in.interval.df$GENENAME)) > 0){
		             write.table(data.frame(CHR=chr, 
				                    CNVTYPE=intervals.df$TYPE[i],
					           CNVSTARTPOS=cnv.start.pos,
						  CNVENDPOS=cnv.end.pos, 
						 GENENAME=unique(genes.in.interval.df$GENENAME)), 
					                 file=genes.in.cnvs.file, col.names=F, row.names=F, 
							                 append=T, quote=F, sep="\t")
		        }
		      }
}

Test <- function(){
    #intervals.file <-"/ifs/scratch/c2b2/ngs_lab/db2175/Projects/PersonalizedCancer/Data/GT/regions.final.txt"
    intervals.file <- "/Users/el2622/Documents/cancerdata/cnvintervalfile.txt"
  refseq.file <- "/Users/el2622/Documents/cancerdata/TruSeq_exome_targeted_regions.hg19.bed.chr"
    #genes.in.cnvs.file <- "/ifs/scratch/c2b2/ngs_lab/db2175/Projects/PersonalizedCancer/Data/GT/regions.final.genes.txt"
    genes.in.cnvs.file <- "/Users/el2622/Documents/cancerdata/geneannotations.txt"
    FindGenesInCNVs(intervals.file, refseq.file,
		                           genes.in.cnvs.file)
}


 args <- commandArgs(TRUE)
 if (is.null(args)){
      print("Please provide parameters")
    exit
     }else{
          print(args)
     }

 intervals.file <- args[1]
 ref.seq.file <- args[2]
 print(paste(args[1],args[2]))
 genes.in.cnvs.file <- paste(intervals.file, "genes.txt",sep="")
 FindGenesInCNVs(intervals.file, ref.seq.file, genes.in.cnvs.file)

