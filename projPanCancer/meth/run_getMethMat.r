P##!/usr/bin/Rscript
#Author: Jing He
#Date:Jul. 17th, 2013 
#Last Updated: Jul. 29rd, 2013
#Usage: 
#Description: Input: tumor acronym, prepare all compressed files and unzip them; output: a txt file with genes and Fisher's Method chi-square statistics
#TODO: check the directory structure, make sure all sub direct have the same structure.

##----------------------------
#receive tumor acrynom from command
# wd <- "/Volumes/ac_lab/jh3283/SCRATCH/projPanCancer/"

args <- commandArgs(TRUE)
if (is.null(args)){
  print("Please provide parameters")
  exit
}else{
  print(args)
}

rootd <- "/Volumes/ac_lab/jh3283/"
tumAcro <- "prad"
# tumAcro <- args[1]

sd <- paste(rootd,"SCRATCH/projPanCancer/meth/",sep="")
cwd <- paste(sd,"/",tumAcro,sep="")
setwd(cwd)
#system("for f in $(ls jhu-usc.edu_COAD.HumanMethylation*tar.gz) ; do tar -xvf ${f} -C . ; done")


# source("http://bioconductor.org/biocLite.R")
# biocLite("GenomicFeatures")
require(GenomicFeatures)
print("loading expression matrix...")
fnexp <- paste("~/SCRATCH/projPanCancer/PANCANCER/",tumAcro,"/",tumAcro,"_tcga_rnaseq.exp",sep="")
dset <- read.table(fnexp, header=T)

# source("/ifs/scratch/c2b2/ac_lab/jh3283/projPanCancer/meth/getMethMat.r")
source(paste(rootd,"scripts/projPanCancer/meth/getMethMat.r",sep=""))
##----------------------------
#For the older file structure two layer

subd <- gsub(".tar.gz","",list.files(pattern="tar.gz"))
for (i in 1:length(subd))
{	
	subdf <- paste(cwd,"/",subd[i],sep="")
	setwd(file.path(cwd, subd[i]))
	subdfs <- dir(path=subdf, pattern=toupper(tumAcro))
	dpattern <- "HumanMethylation"
	barcode=NULL
	print(paste("#Total number of folders:",length(subd)))	
	if(i == 1) {
		meth <- methRead(dpattern=dpattern,barcode=NULL)$meth
		# linkage <- methRead(dpattern=dpattern,barcode=NULL)$linkage
		print (paste("working folder:",subd[i],"number of samples:",length(dir(pattern=dpattern)),sep="\t"))
	}else{
	print (paste("working folder:",subd[i],"number of samples:",length(dir(pattern=dpattern)),sep="\t"))
	tmp <- methRead(dpattern=dpattern,barcode=NULL)$meth
	genes <- unique(c(rownames(meth), rownames(tmp)))
	meth <- cbind(meth[match(genes, rownames(meth)),],tmp[match(genes, rownames(meth)),])
	}
	if(i %% 5 ==0 ) save(meth,  file=paste(cwd,"/",tumAcro, "-meth.rda", sep=""))
	# linkage <- methRead(dpattern=dpattern,barcode=NULL)$linkage
}
save(meth,  file=paste(ofile, "-meth.rda", sep=""))
print("-------------The End----------------")

##----------------------------
#for new directory structure!!
setwd(paste(cwd,"/level_3",sep=""))
dpattern <- "HumanMethylation"
barcode=NULL
meth <- methRead(dpattern=dpattern,barcode=NULL)$meth
save(meth,  file=paste(ofile, "-meth.rda", sep=""))
print("-------------The End----------------")
