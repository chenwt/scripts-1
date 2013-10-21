e!/usr/bin/Rscript
#Author: Jing He
#Date: 
#Last Updated:
#Usage:
#Description:


##--------------------------------------
setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projPanCancer/myethy/")
# prepare data
cwd <- "/ifs/scratch/c2b2/ac_lab/jh3283/projPanCancer/myethy/"
tumAcro <- "brca"
fnexp <- paste("~/SCRATCH/projPanCancer/PANCANCER/",tumAcro,"/",tumAcro,"_tcga_rnaseq.exp",sep="")
exp <- read.table(fnexp)
fnmethy <- paste(cwd,tumAcro,"/",tumAcro,"_genelist_integrate_15_FisherM.txt",sep="")
methyg <- read.table(fnmethy,header=T,sep="\t")
row.names(methyg) <- methyg[,1]

