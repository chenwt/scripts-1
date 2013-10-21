##!/usr/bin/Rscript
#Author: Jing He
#Date:Sep 19,2013 
#Last Updated:
#Usage:
#Description: TARGET AML methylation data

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("lumi")
# source("http://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation27k.db")
# biocLite("IlluminaHumanMethylation27kmanifest")

require(lumi)

rootd = "/ifs/home/c2b2/ac_lab/jh3283/"
cwd = paste(rootd,"~/SCRATCH/projAML/methyl/",sep="")


datad = paste(cwd,"data/",sep="")
outputd = paste(cwd,"result/",sep="")
tempd = paste(cwd,"temp/",sep="")
scriptd =paste(rootd,"~/scripts/projAML/methyl/",sep="")


setwd(cwd)

fns = sort(list.files(path="data/"))


##----------------------------
#test
fnD = fns[grep("Dx",fns)[1]]
fnR = fns[grep("Rm",fns)[1]]
setwd(tempd)
fnR <-"PAEABM_Rm_trimed.txt"
fnD <- "PAEABM_Dx_trimed.txt"
methD = lumiMethyR(paste("temp/",fnD,sep=""), lib="IlluminaHumanMethylation27k.db")
methR = lumiMethyR(paste("temp/",fnR,sep=""), lib="IlluminaHumanMethylation27k.db")
