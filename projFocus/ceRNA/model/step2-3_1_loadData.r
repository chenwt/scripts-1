setRootd = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysnameß']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/projFocus/"
  }
  return(rootd)
}
rootd = setRootd()
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

wd          = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/02022014/model/CHEK1-temp/"
require(gplots)
##---init
gene = tail(unlist(strsplit(wd,split="/|-",perl=T)),2)[1]
setwd(wd)
gene = 'KIF23'
fsmps = "samples.txt"
fexp = "exp.mat"
fsom = "som.mat"
fcnv = "cnv.mat"
fsnp = "snp.mat"
fmeth = "meth.mat"

trimSmpName = function (c) {
  return(vapply(c,FUN=function(x){substr(x,6,15)},'a'))
}

##loading data
rawExp = read.delim2(fexp)
smpExp = colnames(rawExp)
colnames(rawExp) = smpExp
gExp = rownames(rawExp)

rawSom = read.delim2(fsom)
smpSom = trimSmpName(colnames(rawSom))
colnames(rawSom) = smpSom
gMut = rawSom[,1:4]

rawCnv = read.delim2(fcnv)
smpCnv = trimSmpName(colnames(rawCnv))
colnames(rawCnv) =  smpCnv
gCnv = as.character(unlist(rawCnv[,1]))

rawSnp = read.delim2(fsnp)
smpSnp =  trimSmpName(colnames(rawSnp))
colnames(rawSnp) = smpSnp
gSnp = rawSnp[,1:2]


rawSmps = vapply(unlist(read.table("samples.txt")),as.character,'a')
smps = Reduce(intersect,list(rawSmps,smpsCnv,smpSom,smpSnp,smpExp))

#----
dataExp = apply(subset(rawExp,select=smps),2,as.numeric  )
dataExp = t(dataExp[which(rowSums(dataExp)>0),])
dataExp = dataExp[apply(dataExp!=0,1,any),,drop=FALSE]
colnames(dataExp) = gExp 


idx = (grep(gene,gMut[,1],perl=T))
dataSom = apply(subset(rawSom[idx,],select=smps),2,as.numeric  )
dataSom = t(dataSom[apply(dataSom!=0,1,any),,drop=FALSE])
colnames(dataSom) = gSom

dataCnv = subset(rawCnv,select=smps)
dataCnv = t(dataCnv[apply(dataCnv!=0,1,any),,drop=FALSE])
colnames(dataCnv) = gCnv 

idx = (grep(gene,gSnp[,1],perl=T))
dataSnp = apply(subset(rawSnp[idx,],select=smps),2,as.numeric)
dataSnp = t(dataSnp[apply(dataSnp!=0,1,any),,drop=FALSE])
colnames(dataSnp) = gSnp 
