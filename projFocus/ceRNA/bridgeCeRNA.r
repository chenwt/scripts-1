##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: gene list > 
#output: <network file> <network plot>
#Description: this file was created for projFocus, ceRNA, used after running grpLassoSNP.r to map genes to interactom
#TODO: 

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/interactom/")
  rootd = "/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/interactom/")
  rootd = "/ifs/scratch/c2b2/ac_lab/jh3283/"
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
# setwd(system("pwd",intern=T))
# cwd         = getwd()
# filepref    = args['file']
# cutoff      = as.numeric(args['cut'])
# output      = paste(args['out'],".rda",sep="")
# print(paste("current working directory:",cwd))

#-------------functions---
jxy = function(x,y){return(paste(x,y,""))}
#---------------
cwd = paste(rootd,"projFocus/ceRNA/result/interactom/",sep="")
require(igraph)
library(help="igraph")

data = read.delim(paste(cwd,"RANBP9.ceRNAnet.txt",sep=""))
dataNet = graph.data.frame(data,directed=T)
plot(dataNet,edge.color="gray")
