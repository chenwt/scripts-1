#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: CG_cnv.mat CG_methdiff.mat CG_som.mat 
#output: <file: gene Sample list, sample size > 10, all gemonic intact cancer genes>
#Description: this file was created for projFocus, ceRNA, used in step1 to eliminate tumor samples
#     this version includes somatic mutation matrix, next step is to do DEG for each gene on final sample

# ###-----local test
# cnv         = jxy(rootd, "result/02022014/cnv/brca_cnvTumor_level2_combinedCG_02242014.mat") 
# meth        = jxy(rootd, "result/02022014/meth/brca_methTumor_combinedCG_03032014.mat_diffMeth.mat")
# som	        = jxy(rootd, "/result/02022014/som/brca_somTumor_combinedCG_20140301.mat.mat")
# outfile	    = jxy(rootd, "/result/02022014/geneSamples/brca_gslist_combinedCG_CnvMethSomFree_2014-03-03.txt")
# ###------local test data

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
  rootd = "//ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  figd = "/ifs/data/c2b2/ac_lab/jh3283/projFocus/report/topDown_02042014/fig/"
}


#-----funcsStart
barcode2pid = function(x){ return(substr(x,6,15)) }

#----funcEnd
 
args = getArgs()
usage = "Usage: step1-6_getCNVMethFreeSample.r --cnv <cnv.mat> --meth  <methDiff.mat> --som  <somCGtarget.mat> --out <outGeneSampleCNVMethFree.txt >"
example = "Usage: step1-6_getCNVMethFreeSample.r --cnv <cnv.mat> --meth  <methDiff.mat> --som <somCGtarget.mat> --out <outGeneSampleCNVMethFree.txt >"
if(length(args) < 4 || is.null(args)){
   print(usage)
   print(example)
   print(args)
   stop("Input parameter error!")
}else{
   print(args)
}

sampleCnt_cut = 10
setwd(system("pwd",intern=T))
cwd         = getwd()
cnv	        = args['cnv'] 
meth	      = args['meth']
som	        = args['som']
outfile	    = args['out'] 
print(paste("current working directory:",cwd))

##--load data
dcnv = read.delim2(cnv)
dcnv = dcnv[,-ncol(dcnv)]
rownames(dcnv) = dcnv$barcode
dcnv = apply(dcnv[,-1],c(1,2),function(x){ifelse(as.numeric(x)>0,1,0)})
colnames(dcnv) = vapply(colnames(dcnv),barcode2pid,'a')

dmeth = read.delim2(meth,header=T)
dmeth = apply(dmeth[, -ncol(dmeth)],c(1,2),as.numeric)
colnames(dmeth) = vapply(colnames(dmeth),barcode2pid,'a')

dsom = read.delim(som, header=T)
# dupSmp = c("TCGA.B6.A0RE.01A.11D.A060","TCGA.BH.A0H6.01A.21D.A060")
# dsom = dsom[,-sapply(dupSmp, FUN=function(x){grep(x,colnames(dsom))})]
rownames(dsom)= dsom[,1]
dsom = apply(dsom[, -c(1,ncol(dsom))],c(1,2), as.numeric)
smpsom = vapply(colnames(dsom),barcode2pid,'a')
colnames(dsom) = smpsom
# hist(colSums(dsom),freq=T)

smpmeth = colnames(dmeth)
smpcnv = colnames(dcnv)
smpsom = colnames(dsom)
commSmps = intersect( smpmeth, intersect( smpcnv, smpsom) )
dmeth = subset(dmeth,select=commSmps)
dcnv  = subset(dcnv,select=commSmps)
dsom = subset(dsom,select=commSmps)

gmeth = rownames(dmeth)
gcnv = rownames(dcnv)
gsom = rownames(dsom)
commg = intersect( gmeth, intersect( gcnv, gsom) )
dmeth = dmeth[ gmeth %in% commg, ] 
dcnv  = dcnv[ gcnv %in% commg, ] 
dsom  = dsom[ gsom %in% commg, ] 
# ##compare.list(list(rownames(dcnvtemp)), list(rownames(dmethtemp)))
comm = dcnv + dmeth + dsom

pdf(paste(rootd,"/report/topDown_02042014/fig/combinedCancerGeneCnvMethSom",Sys.Date(),".pdf",sep=""),height=10)
require(gplots)
heatmap.2(dcnv,dendrogram='none',col=c("white","blue"),trace='none',labRow="")
heatmap.2(dmeth,dendrogram='none',col=c("white","green"),trace='none',labRow="")
heatmap.2(dsom,dendrogram='none',col=c("white","red"),trace='none',labRow="")
heatmap.2(ifelse(comm==0,0,1),dendrogram='none',col=c("white","black"),trace='none',cexRow=0.5)
dev.off()

result = matrix(0,ncol=length(commSmps),nrow=length(commg))
header = paste("gene","samples_CnvMethSomFree",sep="\t")
write(header,outfile)
for (i in 1:nrow(comm)){
    if ( length(commSmps[which(comm[i,]==0)]) >= sampleCnt_cut ) {
        line = paste(commg[i],paste(commSmps[which(comm[i,]==0)], collapse = ';'),sep="\t")
        write(line,outfile,append=T)
    }
}

