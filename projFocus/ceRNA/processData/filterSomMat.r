#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#TODO: underdevelopmetn

barcode2pid = function(x){ return(substr(x,6,15)) }
jxy = function(...){
  ss = unlist(list(...))
  temp = ""
  res = paste(ss, collapse="")
  return(res)
}

sysInfo = Sys.info()
if(sysInfo['sysname']=="Darwin" ){
  rootd = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/"
  source("/Volumes/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/")
}
else if(sysInfo['sysname']=="Linux" ){
  setwd("/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/meth")
  rootd = "/ifs/ifs/data/c2b2/ac_lab/jh3283/projFocus/"
  source("/ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/projFocusCernaFunctions.R")
  print("working from Linux")
}


s3 = jxy()
somlevel3 = read.delim2("")