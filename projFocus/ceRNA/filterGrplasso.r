##!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: <file: output files prefix from grpLassoSNP.r > 
#output: <r data file: >
#Usage: Rscript grplassoSNP.r input
#Description: this file was created for projFocus, ceRNA, used after running grpLassoSNP.r to filter significant gene-snps 
#TODO: develop plot option!

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test")
}else if(sysInfo['sysname']=="Linux" ){
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
  setwd("/ifs/scratch/c2b2/ac_lab/jh3283/projFocus/ceRNA/result/grpreg/test")
}

###-----------------functions-----------------

getGeneCoeff = function(file,cut){
  temp = readLines(file,n=1)
  temp = unlist(strsplit(unlist(strsplit(temp,,split="\t"))[2],":"))
  if (as.numeric(tail(temp,1)) < cut){
    data = read.delim(file,header=F,skip=1,row.names=1)
    rs = row.names(data)
    data = unlist(data)
    names(data) = rs
    res = data
#      res = list(data)
#     names(res) = temp[1]
    return(res)
  }else{
    res=NULL
    names(res) = NULL
    return(res)
  }
}
getGeneVarNet = function(filelist,cut){
  cntGene = length(filelist)
  res = (rep(NA,cntGene))
  res = sapply(filelist,FUN=function(x){getGeneCoeff(x,cut)})
  #this line assume there is no dot in a gene name,the inputfile format
  names(res) = vapply(names(res),FUN=function(x){ gsub(".txt","",tail(unlist(strsplit(x,"_")),1))},'a') 
  res[sapply(res, is.null)] <- NULL
  return(res)
}

###------------------functions---end---------


args = getArgs()
usage = "Usage: Rscript filterGrplasso.r --file <grplasso coeff prefix> --cut <0.05, smallest is 1/nperm, default 0.05> --out <filtered rda file name> "
example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/filterGrplasso.r --file grplasso_coeff --cut 0.05 --out gcGenes_GeneVarNet"
if(length(args) < 3 || is.null(args)){
  print(usage)
  print(example)
  print(args)
  stop("Input parameter error!")
}else{
  print(args)
}

setwd(system("pwd",intern=T))
cwd         = getwd()
filepref    = args['file']
cutoff      = as.numeric(args['cut'])
output      = paste(args['out'],".rda",sep="")
print(paste("current working directory:",cwd))

##-----------------------test-##----------------------------
# setwd(system("pwd",intern=T))
# cwd         = getwd()
# filepref    = "grplasso_coeff_grplasso_"
# cutoff      = 0.05
# output      = "gcGrpVarGenNet.rda"
# print(paste("current working directory:",cwd))

##----------------------------##----------------------------
fileList    = dir(pattern=filepref)
test = getGeneCoeff(file=fileList[2],cut=0.05 )
geneVarNet  = getGeneVarNet(fileList,cutoff)
save(geneVarNet,file=output)
write.table(as.matrix(names(geneVarNet)),paste(args['out'],".genelist.txt",sep=""),
            col.names=F,row.names=F,quote=F)
print("#---DONE---")
