#!/bin/Rscirpts
#cat AC1SNPs.txt | awk '{print $1"\t"$2+1"\t"$3}' > AC1SNPs_3.txt 
#example: Rscript ~/scripts/projNET/run_PlotMA.r 13 19

args <- commandArgs(TRUE)

if (is.null(args)){
	print("Please provide parameters: ac#,and patient snp file name")
	exit
}else{
	t.num <- args[1]
	n.num <- args[2]
	print(c(t.num,n.num))
}

SCRIPTSDIR <- Sys.getenv("SCRIPTS")
###------------------------------ preparing filenames
NETDIR <- Sys.getenv("NETDIR")
if(is.null(NETDIR)) {exit}
cwd <- paste(NETDIR,"AC",t.num,"/res/",sep="")
tumor.file <- paste(NETDIR,"AC",t.num, "/ac",t.num,"_fpkm_gene",sep="") 
normal.file <- paste(NETDIR,"AC",n.num, "/ac",n.num,"_fpkm_gene",sep="") 
del_reg.file <- paste(NETDIR, "AC3somatic_corrected_NoModLong-3-type.seg",sep="")

print(c("Working directory: ",cwd))
print(c("Scripts directory: ", paste(SCRIPTSDIR,"projNET/MAplot.r",sep="")))
print(c("Tumor fpkm file:", tumor.file))
print(c("Normal fpkm file:", normal.file))
print(c("Del_region file:", del_reg.file))

###------------------------------ loading packages and source code, run functions
setwd(cwd)
source(paste(SCRIPTSDIR,"projNET/PlotMA.r",sep=""))
require(ggplot2)

MAplot(tumor.file,normal.file,t.num,n.num,del_reg.file)








