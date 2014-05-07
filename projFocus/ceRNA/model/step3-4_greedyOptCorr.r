# CWD = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/"
# setwd(CWD)

##test files------
# expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/CEP55_exp"
# mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/CEP55_regMut"
# figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/May2014/fig/"
#  output = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_sigMut.txt"
##test end------
usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,3,6)],collapse="-")

setRootd  = function(){
  sysInfo = Sys.info()
  if(sysInfo['sysname']=="Darwin" ){
    print("working from MacOS")
    rootd = "/Volumes/ifs/home/c2b2/ac_lab/jh3283/"
  }else if(sysInfo['sysname']=="Linux" ){
    print("working from Linux")
    rootd = "/ifs/home/c2b2/ac_lab/jh3283/"
  }
  return(rootd)
}
rootd     = setRootd()
source(paste(rootd,"/scripts/projFocus/ceRNA/projFocusCernaFunctions.R",sep=""))
source(paste(rootd,"/scripts/myR/jingGraphic.R",sep=""))

# #----getting command line parameters
args = getArgs()
if(length(args) < 1 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop(paste(error,"wrong input parameter!"))
}

setwd(system("pwd",intern=T))
expfile     = args$exp
mutfile     = args$mut
output      = args$output
figd        = paste(rootd, "/DATA/projFocus/report/May2014/fig/", sep="")
print("###------------------------")
print(paste("inputfile",expfile, mutfile))
# print(paste("outputfile",output))
#---init

###---func
fisherZ = function(r) {
  if (abs(r) == 1){
    return(r * 4.95)
  }else{
    return(1/2*log((1 + r )/(1-r)))
  }
}

z2corr = function(z) return((exp(2 * z) -1)/(exp(2 * z) +1))

corrDiff = function(z1,z2,n1,n2)  {
  if (n1 >=5 & n2 >=5){
  z = (z2 - z1) / sqrt( 1/(n2-3) + 1/(n1-3) ) 
  return(2*(1-pnorm(abs(z))))
  }else{
    return(NA)
  }
}

calCorr = function(expD, mutD, tarExp){
  sumExp = colSums(expD * mutD)
  sumExp = sumExp[which(sumExp != 0)]
  numCol = length(sumExp)
  ##need further optimization  
  resCorr = fisherZ(cor(sumExp, tarExpD[names(sumExp)])) 
  return(list(zs=resCorr, n = numCol))
}

corrOpt_binflip = function(mutD, regExpD, tarExpD, corr_init = corr_full, tol = 0){
  mutInd = which(mutD>0,arr.ind=T)
  numSmp = length(unique(mutInd[,2]))
  numReg = length(unique(mutInd[,1]))
  if (numSmp <= 3 | numReg <=3){
    return(list(corr_opt =NA, sample =unique(mutInd[,2]) , cerna=unique(mutInd[,1]), corrDiff_pval=NA, mutD=mutD[rowSums(mutD)!=0,(colSums(mutD)!=0)]))
  }else{
    numMut = nrow(mutInd)
    mutInd = mutInd[sample(numMut),]
    mut_temp = mutD
    corr_prev = corr_init
    
  
    for (i in 1:numMut){
      id_flip = mutInd[i,]
      mut_temp[id_flip[1], id_flip[2]] <- 0
      corr_temp = calCorr(regExpD,mut_temp,tarExpD)
      if (is.na(corr_temp$zs)) {break}
      if ( abs(corr_temp$zs) - abs(corr_prev$zs) < tol){
        mut_temp[id_flip[1], id_flip[2]] <- 1
        corr_temp = list(corr=0,n=0)
      }else{
        corr_prev = corr_temp
        corr_temp = list(corr=0,n=0)
      }
    }
    corrDiffPval = corrDiff(corr_prev$zs, corr_full$zs,corr_prev$n,corr_full$n)
    
    smpOpt = colSums(mut_temp)[colSums(mut_temp)>0]
    numSmpOpt = length(smpOpt)
    regOpt = rowSums(mut_temp)[rowSums(mut_temp)>0]
    numRegOpt = length(regOpt)
    corr_opt = corr_prev
    corr_opt$corr <- (exp(2 * corr_prev$zs) -1)/(exp(2 * corr_prev$zs) +1)
    return(list(corr_opt =corr_opt, sample = smpOpt, cerna=regOpt,corrDiff_pval=corrDiffPval, mutD=mut_temp))
  }
}

permuMutD = function(regmutD){
  colname_orig = colnames(regmutD)
  rowname_orig = row.names(regmutD)
  mutD_perm = regmutD
  mutD_perm = mutD_perm[sample(nrow(mutD_perm)),]
  mutD_perm = mutD_perm[,sample(ncol(mutD_perm))]
  colnames(mutD_perm) <- colname_orig
  rownames(mutD_perm) <- rowname_orig 
  return(mutD_perm)
}

doPermu = function(regExpD, regmutD, tarExpD, nperm = 1000){
  corr_perm = rep(0, nperm)
  for (i in 1:nperm){
    corr_perm[i]  = calCorr(expD=regExpD,mutD=permuMutD(regmutD),tarExp=tarExpD)$zs
  }
  return(corr_perm)
}

plot_perm = function(corr_perm, resCorr, corr_full){
  hist(corr_perm, col="lightgray", xlim=c(-0.8,0.8),border="gray",
       xlab = "permut correlation", main = " correlation from permutation")
  abline(v=resCorr$corr_opt$zs,col="red",lwd=4)
  abline(v=corr_full$zs,col="orange",lwd=4)
  abline(v=corr_total, col="green", lwd = 4,lty=2)
}

removeZeor = function(arr) return(arr[rowSums(arr)!=0, colSums(arr)!=0])

writeOut = function(file, mutD ){
  outHeader = paste(c("tgene","totalSmp","mutSmp","selectMutSmp","totalReg","mutReg","selectMutReg","totalCorr","fullCorr","optCorr","permuCorr","selectGenes","selectSamples"),collapse="\t")
  resMut = removeZeor(removeZeor(mutD))
  outGene = rownames(resMut)
  outSmp = colnames(resMut)
  outSmpMap = vector(length=length(outGene))
  for (i in 1:nrow(resMut)){
    outSmpMap[i] = paste("[",paste(vapply(which(resMut[i,]!=0, arr.ind=T),function(x){outSmp[x]},'a'),collapse=""), "]",sep="")
  } 
  outRecord = paste(c(tgene, ncol(expD), ncol(resMut), NA, nrow(mutD), nrow(resMut),NA, NA, NA,NA, NA, 
                      paste(outGene, collapse=";"),
                      paste(outSmpMap,collapse=";") ), collapse="\t")
  
  write.table(outHeader, file=file, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  write.table(outRecord, file=file, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
}
###----func

expD = read.table(expfile,sep="\t",header = T)
rownames(expD) = expD[,1]
expD = apply(expD[,-1],c(1,2),as.numeric)
exp_min = min(apply(expD,1,min))
expD = expD - exp_min

mutD = read.table(mutfile,sep="\t", header = T)
rownames(mutD) = mutD[,1]
mutD = apply(mutD[,-1], c(1,2), as.numeric)

tarExpD = expD[1,]
tgene = rownames(expD)[1]
tarExpD = tarExpD[order(tarExpD)]

### mutation sample, regulator number check
if (is.null(removeZeor(mutD))){
  flag = FALSE
}else if (NROW(removeZeor(mutD)) < 2){
  flag = FALSE
}else if (NCOL(removeZeor(mutD)) < 5) {
  flag = FALSE
} else if (nrow(expD) < 3 ) {
  flag = FALSE  
}else {
  flag = TRUE
}

if (flag) {
  writeOut(output, mutD)
  write.table("#no enough samples/regulators", file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  q(save = "no")
}

regExpD = subset(expD[-1,],select= names(tarExpD))
sumExpD = colSums(regExpD)
corr_total = cor(sumExpD,tarExpD)
sumExpD = colSums(tarExpD * mutD)

corr_full = calCorr(regExpD,mutD,tarExpD)

if (NCOL(removeZeor(mutD)) < 5 ) {
  outHeader = paste(c("tgene","totalSmp","mutSmp","selectMutSmp","totalReg","mutReg","selectMutReg","totalCorr","fullCorr","optCorr","permuCorr","selectGenes","selectSamples"),collapse="\t")
  resMut = removeZeor(removeZeor(mutD))
  outGene = rownames(resMut)
  outSmp = colnames(resMut)
  outSmpMap = vector(length=length(outGene))
  for (i in 1:nrow(resMut)){
    outSmpMap[i] = paste("[",paste(vapply(which(resMut[i,]!=0, arr.ind=T),function(x){outSmp[x]},'a'),collapse=""), "]",sep="")
  } 
  outRecord = paste(c(tgene, ncol(expD), ncol(resMut), NA, nrow(mutD), nrow(resMut),NA, corr_total, corr_full,NA, NA, paste(outGene, collapse=";"), paste(outSmpMap,collapse=";")), collapse="\t")
  write.table(outHeader, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  write.table("#no enough mutated samples", file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  write.table(outRecord, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  q(save="no")
}

nperm = 1000
tolVec = seq(-0.005,0.005,by=0.001)
nTol = length(tolVec)

##select tolence
iterTolSum = as.data.frame(matrix(NA, nrow=nTol,ncol=5))
colnames(iterTolSum) = c("tol", "optcorr", "optSmpn","pvalPerm", "pvalfull")
for (i in 1:nTol ){
  tol = tolVec[i]
  resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD,corr_full, tol = tol)
  corr_perm = doPermu(regExpD, resCorrOpt$mutD, tarExpD)
  pval_perm = max(length(corr_perm[corr_perm > resCorrOpt$corr_opt$zs])/nperm, 0.00001)
  pval_full = resCorrOpt$corrDiff_pval
  iterTolSum[i, ]  = c(tol = tol, optcorr=z2corr(resCorrOpt$corr_opt$zs), optSmpn = round(resCorrOpt$corr_opt$n), pvalPerm = pval_perm, pvalfull=pval_full)
}

out = paste(output, ".iterTols.temp", sep="")
write.table(iterTolSum, out, col.names=F, sep="\t",quote=F)
pvalCut = 0.05;
tol = iterTolSum[do.call(order, list(iterTolSum[,5],iterTolSum[,4],iterTolSum[,3])),]
tol = na.omit(tol)
if (length(tol[,1]) > 0 ){
  if(tol[1,4] <= pvalCut & tol[1,5] <= pvalCut){
    tol = tol[which.max(tol[,2]),1]
  }else{
    outHeader = paste(c("tgene","totalSmp","mutSmp","selectMutSmp","totalReg","mutReg","selectMutReg","totalCorr","fullCorr","optCorr","permuCorr","selectGenes","selectSamples"),collapse="\t")
    outRecord = paste(c(tgene, ncol(expD), corr_full$n, iterTolSum[1,3], nrow(expD) - 1, nrow(removeZeor(mutD)), NA, corr_total, tol[1,4], iterTolSum[1,2], NA,NA, NA), collapse="\t")
    write.table(outHeader, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
    write.table("#optimization pvalue not significant", file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
    write.table(outRecord, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
    q(save="no")
  }
}else{
  outHeader = paste(c("tgene","totalSmp","mutSmp","selectMutSmp","totalReg","mutReg","selectMutReg","totalCorr","fullCorr","optCorr","permuCorr","selectGenes","selectSamples"),collapse="\t")
  outRecord = paste(c(tgene, ncol(expD), corr_full$n, NA, nrow(expD) -1 , nrow(removeZeor(mutD)), NA, corr_total,tol[1,4], NA, NA,NA, NA), collapse="\t")
  write.table(outHeader, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  write.table("#no enough mutated samples", file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  write.table(outRecord, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
  q(save="no")
}

##random initiation
nIter = 100
resIterMutD = matrix(0, nrow=nrow(resCorrOpt$mutD),ncol = ncol(resCorrOpt$mutD))
iterRansInitSum = as.data.frame(matrix(NA, nrow=nTol,ncol=5))
colnames(iterRansInitSum) = c("#inter", "optcorr", "optSmpn","pvalPerm", "pvalfull")
pvalCut = 0.05; corr_max = 0.0
for (i in 1:nIter ){
  resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD,corr_full, tol = tol)
  corr_perm = doPermu(regExpD, resCorrOpt$mutD, tarExpD)
  pval_perm = max(length(corr_perm[corr_perm > resCorrOpt$corr_opt$zs])/nperm, 0.00001)
  pval_full = resCorrOpt$corrDiff_pval
  if(!is.na(pval_full) & !is.na(pval_perm) & pval_full > pvalCut & pval_perm > pvalCut & z2corr(resCorrOpt$corr_opt$zs) > corr_max){
    corr_max =  z2corr(resCorrOpt$corr_opt$zs)
    resMutOptD= resCorrOpt$muD
  }
    resIterMutD = resIterMutD + resCorrOpt$mutD
    iterRansInitSum[i, ]  = c(i, optcorr=z2corr(resCorrOpt$corr_opt$zs), optSmpn = round(resCorrOpt$corr_opt$n), pvalPerm = pval_perm, pvalfull=pval_full)

}

out = paste(output, ".randomInit.temp", sep="")
write.table(iterRansInitSum, out, row.names=F,sep="\t",quote=F)

ncut = 30
resMut = apply(resIterMutD, c(1,2), function(x){ifelse(x > ncut, 1, 0)})
corr_final = calCorr(expD=regExpD, mutD=resMut, tarExp=tarExpD)
resMut = removeZeor(resMut)
outGene = rownames(resMut)
outSmp = colnames(resMut)
outSmpMap = vector(length=length(outGene))
for (i in 1:nrow(resMut)){
  outSmpMap[i] = paste("[",paste(vapply(which(resMut[i,]!=0, arr.ind=T),function(x){outSmp[x]},'a'),collapse=""), "]",sep="")
}

outNumbers = c(tgene, ncol(expD), corr_full$n, corr_final$n, nrow(expD) -1 , nrow(removeZeor(mutD)), nrow(resMut), corr_total, z2corr(corr_full$zs), z2corr(corr_final$zs), z2corr(mean(corr_perm)))
outSmpMap = paste(outSmpMap,collapse=";")
outGene = paste(outGene, collapse=";")
outRecord = paste(c(outNumbers,outGene,outSmpMap),collapse="\t")
outHeader = paste(c("tgene","totalSmp","mutSmp","selectMutSmp","totalReg","mutReg","selectMutReg","totalCorr","fullCorr","optCorr","permuCorr","selectGenes","selectSamples"),collapse="\t")
write.table(outHeader, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
write.table(outRecord, file=output, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)

print("[END]")
