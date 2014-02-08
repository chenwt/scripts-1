##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file1: tumor methlation matrix> <file2: normal methylation matrix> 
#output: <file: tumor sample with methylation relative to population mean>
#Description: this file was created for projFocus, ceRNA, used in step1 to eliminate tumor samples
#TODO: 

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/")
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  figd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth")
  rootd = "/ifs/scratch/c2b2/ac_lab/jh3283/"
  figd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}
# 
# args = getArgs()
# usage = "Usage: Rscript bridegCeRAN.r --file <gene.list>  "
# example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterGrplasso.r --file grplasso_coeff --cut 0.05 --out gcGenes_GeneVarNet"
# if(length(args) < 3 || is.null(args)){
#   print(usage)
#   print(example)
#   print(args)
#   stop("Input parameter error!")
# }else{
#   print(args)
# }
# 

cnv = "cnv/brca_cnvTumor_level2_ucCeRNETCancerGene_02062014.mat"
meth = "meth/brca_methTumor_level3_02072014.mat_diffMeth.mat"
dcnv = read.delim2(cnv)
dcnv = dcnv[,-ncol(dcnv)]
rownames(dcnv) = dcnv$barcode
dcnv = apply(dcnv[,-1],c(1,2),function(x){ifelse(as.numeric(x)>0,1,0)})
temp = vapply(colnames(dcnv),FUN=function(x){substr(x,0,19)},'a')
colnames(dcnv) = temp

dmeth = read.delim2(meth,header=T)
dmeth = apply(dmeth[, -ncol(dmeth)],c(1,2),as.numeric)
temp = vapply(colnames(dmeth),FUN=function(x){substr(x,0,19)},'a') 
colnames(dmeth) = temp
  
sampleComm = colnames(dcnv)[colnames(dcnv) %in% colnames(dmeth)]
geneComm = rownames(dcnv)[rownames(dcnv) %in% rownames(dmeth)]
comm = dcnv[geneComm,sampleComm] + dmeth[geneComm,sampleComm]
sampleNumber = 115 - rowSums(comm)
plot(density(sampleNumber),col="blue")
plot(sampleNumber, col="blue")
sampleUnion = unique(c(colnames(dcnv),colnames(dmeth)))
geneUnion = unique(c(rownames(dcnv),rownames(dmeth)))
result = matrix(0,ncol=length(sampleUnion),nrow=length(geneUnion))

as.matrix(sampleComm)
file = "brca_geneSamplelist_UCceRNET_cancerGene_CNVMethFree_02072014.txt"
header = paste("gene","samples_noCNV_noMeth",sep="\t")
write(header,file)
# for (i in 1:5){
for (i in 1:nrow(comm)){
    line = paste(geneComm[i],paste(sampleComm[which(comm[i,]==0)], collapse = ';'),sep="\t")
    write(line,file,append=T)
}

