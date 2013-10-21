##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep. 26, 2013 
#Last Updated:
#Desp.: calculate the stochastic of one mutation appears in 2 patients
#Input:
#Output:
# TODO: make it reuseable to calculate probability!!
# source("~/scripts/projAML/WXS/get_Pvalue.r")

# rootd = "/ifs/scratch/c2b2/ac_lab/jh3283/"
# cwd = "/Volumes/ac_lab/jh3283/SCRATCH/projAML/WXS/reports/result_labmetgSep24"
cwd ="/ifs/scratch/c2b2/ac_lab/jh3283/projAML/WXS/reports/result_labmetgSep24"
setwd(cwd)
##----------------------------
#load data of raw mutations numbers in each patients
data <- read.table("temp_mutation_count",header=T,fill=T,row.names=1) 
# data <- read.table("temp_mutation_count_filtered",header=T,fill=T,row.names=1) 
colnames(data) <- c("samtools_tn","samtools_rn","gatk_tn","gatk_rn")
colms <- floor(colMeans(data))

##----------------------------
#global variable
# 1. Total #bases in mRNA : 14,881,824,369
# 2. Total #bases in exons : 99,752,470
# 3. Total # bases in RefSeq Genes : 2,011,862,672
ref <- 2011862672
getSample = function(a,b){
	return(sample.int[b,size=a,replace=F])
}

n= 16
getProb = function(n,t) 1 - exp(lfactorial(t) - n * log(t) - lfactorial(t - n))
pvs = sapply(colms,function(x) getProb(n,x))
pvsByP = apply(data,c(1,2),function(x) getProb(n,x))
pvColms <- colMeans(pvsByP)

doPermu = function(colms,np=1000){
	stn <- getSample(colms[1],ref); srn <- getSample(colms[2],ref);
	gtn <- getSample(colms[3],ref); grn <- getSample(colms[4],ref);
	slen <- length(intersect(stn,srn))
	glen <- length(intersect(gtn,grn))

}