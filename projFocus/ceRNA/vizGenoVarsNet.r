##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: gene list > 
#output: <network file> <network plot>
#Description: this file was created for projFocus, ceRNA, used after running grpLassoSNP.r to map genes to interactom
#TODO: 

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test/")
  rootd = "/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/"
  figd = "/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/report/figure/"
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test/")
  rootd = "/ifs/scratch/c2b2/ac_lab/jh3283/"
  figd = "/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/report/figure/"
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
cwd = getwd()
require(igraph)
library(help="igraph")
require(gplots)


####---------------genemoic variancts network---
data = read.delim(paste(cwd,"/grplasso_coeffRANBP9.txt",sep=""),skip=1,header=F)
data = cbind(rep("RANBP9",nrow(data)),data)
dataNet = graph.data.frame(data,directed=F)

V(dataNet)$color = bluered(nrow(data))
V(dataNet)['RANBP9']$color = "darkgreen"
V(dataNet)$size = 5
V(dataNet)['RANBP9']$size = 8
V(dataNet)$label = NA
labs = c("cnv","som","rs16874698", "rs12183724","RANBP9")
V(dataNet)[labs]$label = labs


####----plot
pdf(paste(figd,"/RANBP9_geneticVars_net.pdf",sep=""))
par(mai=c(0,0,1,0))
plot(dataNet,
     layout=layout.fruchterman.reingold, vertex.color=V(dataNet)$color,
     vertex.label.dist=0.2,
     vertex.frame.color="white",
     #      vertex.label.color="black",
     #      vertex.label.font = 2,
     #      vertex.label.cex  = 0.5,
     vertex.label = V(dataNet)$label,
     vertex.size = V(dataNet)$size,
     edge.width = 0.8,
     edge.arrow.size = 0.8,
     edge.arrow.width = 0.2     
)
dev.off()
