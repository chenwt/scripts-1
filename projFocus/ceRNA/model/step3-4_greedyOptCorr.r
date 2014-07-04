#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
rm(list=ls())
usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
typeSelect = "20" #which type of output <all> to select a cutoff <max> to select random optimal(n = 100), or a number to set cutoff value
typeTol = "flex" # < fix> to set tol == 0 , < flex> to select tol adaptively
timeStart = Sys.time()

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
expfile     = args[['exp']]
mutfile     = args[['mut']]
typeTol     = args['ttol']
typeSelect  = args['tsel']
numRandom   = args['nrand']
output      = paste(args['output'], "_", typeTol, "_", typeSelect, "_", numRandom, sep="")
figd        = paste(rootd, "/DATA/projFocus/report/May2014/fig/", sep="")
print("###------------------------")
print(paste("inputfile 1",class(expfile), expfile))
print(paste("inputfile 2",class(mutfile), mutfile))
print(paste("ttol ",class(typeTol), typeTol))
print(paste("tsel",class(typeSelect), typeSelect))
print(paste("output",class(output), output))

print(paste("outputfile",output))
#---init

# ##test files------
# expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp"
# mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp"
# 
# expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/runMay5/flex_max_1000/optCorr.result_SPRY1_exp.temp"
# mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/runMay5/flex_max_1000/optCorr.result_SPRY1_regMut.temp"
# output = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul4/test/test_optCorr_SPRY1.txt"
# 
# figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jul2014/test/"
# output = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_NLK_regMut.temp.txt_test"
# output = paste(output, "_", typeTol, "_", typeSelect, sep="")
# 
# expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/runJul1/test//optCorr.result_FAM126B_exp.temp"
# mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/runJul1/test/optCorr.result_FAM126B_regMut.temp"
# output = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/runJul1/test/test_optCorr_FAM126B.txt"


typeTol = "flex"
typeSelect = "max"
numRandom = 10
# ##test end------


###func-----
fisherZ = function(r) {
  if (is.na(r)){
    return(NA)
  }else if(abs(as.integer(r)) > 1){
    print(paste(ERR, "in correlation:", r))
    return(NA)
  }else if(abs(as.integer(r)) == 1){
    return(r * 4.95)
  }else{
    return(1/2*log((1 + r )/(1-r)))
  }
}

z2corr = function(z) return((exp(2 * z) -1)/(exp(2 * z) +1))

zDiff = function(z1,z2,n1,n2)  {
  if (n1 >=10 & n2 >=10){
    z = (z2 - z1) / sqrt( 1/(n2-3) + 1/(n1-3) ) 
    return(2*(1-pnorm(abs(z))))
  }else if ( ( n1 > 3 & n1 <10) | (n2 >3 & n2 < 10)){
    print( paste(" samples should be larger than 10",n1, n2))
    z = (z2 - z1) / sqrt( 1/(n2-3) + 1/(n1-3) ) 
    return(2*(1-pnorm(abs(z))))
  }else{
    return(NA)
  }
}

calZscore = function(expD, mutD, tarExp){
  if (NROW(expD) != NROW(mutD) | NCOL(expD) != NCOL(mutD)){
    return(list(zs=NA,n=NA))
  }else{
     if (NCOL(expD) == 1)  {
       sumExp = expD * mutD
     }else{
       sumExp = colSums(expD * mutD)
     }
  }
  sumExp = sumExp[which(sumExp != 0)]
  numCol = length(sumExp)
  if (numCol != 0 ) {
    resCorr = fisherZ(cor(sumExp, tarExpD[names(sumExp)])) 
    return(list(zs=resCorr, n = numCol))
  }else{
    return(list(zs=NA, n=NA))
  }
}

removeZeor = function(arr) return(arr[rowSums(arr)!=0, colSums(arr)!=0])

corrOpt_binflip = function(mutD, regExpD, tarExpD, zs_full = zs_full, tol = 0){
  ### one iteration
  ##get mutation indice
  allSmp = colnames(mutD)
  allReg = rownames(mutD)
  
  mutInd = which(mutD>0,arr.ind=T)
  numSmp = length(unique(mutInd[,2]))
  numReg = length(unique(mutInd[,1]))
  
  #QC: check sample size and regulator size 
  if (numSmp < 10 | numReg < 2 ){
    return(list(corr_opt =NA, sample = numSmp , regs = numReg, 
                zDiff_pval=NA, mutD=removeZeor(mutD)))
  }
  
  ##random starting 
  numMut = nrow(mutInd)
  mutInd = mutInd[sample(numMut),]
  mut_temp = mutD
  corr_prev = zs_full
  
  ##flipover forward selection, greedy elemination of mutations
  for (i in 1:numMut){
    id_flip = mutInd[i,]
    mut_temp[id_flip[1], id_flip[2]] <- 0
    corr_temp = calZscore(regExpD,mut_temp,tarExpD)
    if (is.na(corr_temp$zs)) {break}
    if ( abs(corr_temp$zs) - abs(corr_prev$zs) < tol){
      mut_temp[id_flip[1], id_flip[2]] <- 1
      corr_temp = list(corr=0,n=0)
    }else{
      corr_prev = corr_temp
      corr_temp = list(corr=0,n=0)
    }
  }
  
  ## pvalue with full
  zDiffPval = zDiff(corr_prev$zs, zs_full$zs,corr_prev$n,zs_full$n)
  smpOpt = names(colSums(mut_temp)[colSums(mut_temp)>0])
  numSmpOpt = length(smpOpt)
  regOpt = names(rowSums(mut_temp)[rowSums(mut_temp)>0])
  numRegOpt = length(regOpt)
  corr_opt = corr_prev
  return(list(corr_opt = corr_prev, zs_full=zs_full,
              pval_full = zDiffPval,
              sample = smpOpt, regs =regOpt, 
              mutD=mut_temp))

}

permuMutD = function(mutD){
  ### permute row and column labels
  colname_orig = colnames(mutD)
  rowname_orig = row.names(mutD)
  mutD_perm = mutD
  mutD_perm = mutD_perm[sample(nrow(mutD_perm)),]
  mutD_perm = mutD_perm[,sample(ncol(mutD_perm))]
  colnames(mutD_perm) <- colname_orig
  rownames(mutD_perm) <- rowname_orig 
  return(mutD_perm)
}

optCorr = function(mutD, regExpD, tarExpD, nperm) {
  ##simple QC-------
  if (is.null(removeZeor(mutD))){
    flag = FALSE
  }else if (NROW(removeZeor(mutD)) < 2){
    flag = FALSE
  }else if (NCOL(removeZeor(mutD)) < 10) {
    flag = FALSE
  }else {
    flag = TRUE
  } 
  ##----
  regExpD = subset(expD,select = names(tarExpD))[-1,]
  
  if (NCOL(regExpD) == 1){ 
    sumExpD = regExpD 
  }else{
    sumExpD = colSums(regExpD)
  }
  ###------permute actual mutation started he
  
  ### all expression Mreg v.s Etar
  zs_total = list(zs = fisherZ(cor(sumExpD, tarExpD)),
                  n = length(sumExpD) ) 
  
  ### full mutation matrix Mreg * Ereg v.s Etar
  zs_full = calZscore(regExpD,mutD,tarExpD)
  
  if (!flag) {
    #fail qc, exit 
    writeOut(output, mutD, tgene, zs_total, zs_full, tag="failQC")
    q(save = "no")
  }
  
  ##select tolence------
  ### keeping the original method of selecting tol
  nperm = as.integer(numRandom)
  pvalCut = 0.1
  
  if (typeTol == "flex") {
    tolVec = seq(-0.02,0.005,by=0.001)
    nTol = length(tolVec)
    
    iterTolSum = as.data.frame(matrix(NA, nrow=nTol,ncol=5))
    colnames(iterTolSum) = c("tol", "optcorr", "optSmpn",'optRegn', "pvalfull")
    
    for (i in 1:nTol ){
      tol = tolVec[i]
      resCorr_crt = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
      pval_full = resCorr_crt$zDiff_pval
      iterTolSum[i, ]  = c(tol = tol, optcorr=z2corr(resCorr_crt$corr_opt$zs),
                           optSmpn = resCorr_crt$corr_opt$n, 
                           optRegn = length(resCorr_crt$regs),
                           pvalfull = resCorr_crt$pval_full)
    }
    
    
    tol = iterTolSum[which(is.na(iterTolSum$pvalfull)),]
    
    # pvalue QC
    if (NROW(tol)  == 0 ){
      tol = iterTolSum[do.call(order, list(iterTolSum$pvalfull,iterTolSum$optRegn)),][1,1]
      resCorr_crt = corrOpt_binflip(mutD, regExpD, tarExpD, zs_full, tol = tol)
      writeOut(output, mutD, tgene,zs_total=zs_total, zs_full=zs_full, 
               resCorr_crt, tag="failOptTol")
      q(save="no")
    }
    
    tol = tol$tol[ which.max( vapply(tol$pvalfull, FUN=function(x){ifelse(is.na(x), 0, -log(x))},0.1) )[1] ]
    
  }else if (typeTol == "fix") {
    tol = 0
  }else{
    print(paste(ERR, "ttol < flex/fix> :", typeTol))
    q(save="no")
  }
  
  ##random starting-----
  print(paste("Finished Tol time :", format.difftime(difftime(Sys.time(), timeStart)), sep=""))
  
  nIter = as.integer(numRandom)
  pval_full_min <- pval_perm_min <- 1
  
  if ( typeSelect == "max" ) {
    print("seletion using optimal")
    resRandInitOptCorr = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
    
    for (i in 1:nIter){
      resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
      if(!is.na(resCorrOpt$pval_full) & 
           resCorrOpt$pval_full < pval_full_min){
        pval_full_min = resCorrOpt$pval_full
        resRandInitOptCorr = resCorrOpt
      }
    }
    resCorrOpt = resRandInitOptCorr
    resIterMutD = resRandInitOptCorr$mutD
  }else{
    resIterMutD = matrix(0, nrow= nrow(mutD),ncol =ncol(mutD) )
    
    iterRansInitSum = as.data.frame(matrix(NA, nrow=nIter,ncol=5))
    colnames(iterRansInitSum) = c(paste("tol",tol,sep="_"), "optcorr",  
                                  "optSmpn",'optRegn', "pvalfull")
    for (i in 1:nIter){
      resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
      resIterMutD = resIterMutD + resCorrOpt$mutD
      iterRansInitSum[i, ]  = c(iter = paste(tol,i,sep="_"), 
                                optcorr = z2corr(resCorrOpt$corr_opt$zs),
                                optSmpn = resCorrOpt$corr_opt$n, 
                                optRegn = length(resCorrOpt$regs),
                                pvalfull = resCorrOpt$pval_full)
    }  
    if (typeSelect == 'all' | typeSelect == 'max'){
      print("seletion uisng  all ")
      ncut = 0  
      resMut = apply(resIterMutD, c(1,2), function(x){ifelse(x >= ncut, 1, 0)})
      
      cntReg = resCorrOpt$corr_opt$n    
      final_corr = calZscore(expD=regExpD, mutD=resMut, tarExp=tarExpD)
      resCorrOpt <- list( mutD = resMut, corr_opt = final_corr,
                          regs = rownames(removeZeor(resMut)), sample = colnames(removeZeor(resMut)),
                          zs_full = zs_full, pval_full = zDiff(final_corr$zs, zs_full$zs, final_corr$n, zs_full$n))
      
    }else if (length(grep('[0-9.]+', typeSelect)) == 1 & length((grep('[a-z]+', typeSelect))) == 0){
      ncut = as.integer(typeSelect)
      print(paste("seletion using cutoff", ncut))
      
      resMut = apply(resIterMutD, c(1,2), function(x){ifelse(x >= ncut, 1, 0)})
      final_corr = calZscore(expD=regExpD, mutD=resMut, tarExp=tarExpD)
      resCorrOpt <- list( corr_opt = final_corr,
                          regs = rownames(removeZeor(resMut)), sample = colnames(removeZeor(resMut)),
                          zs_full = zs_full, pval_full = zDiff(final_corr$zs, zs_full$zs,final_corr$n,zs_full$n)) 
    }else{
      print( paste(ERR, "tinit should be <max> < all> or <0-100 numbers>", typeSelect))
      #     q(save='no')
    }
  }
  return(corrOptObj=resCorrOpt)
}

plot_perm = function(corr_perm, resCorr, zs_full){
  hist(corr_perm, col="lightgray", xlim=c(-0.8,0.8),border="gray",
       xlab = "permut correlation", main = " correlation from permutation")
  abline(v=resCorr$corr_opt$zs,col="red",lwd=4)
  abline(v=zs_full$zs,col="orange",lwd=4)
  abline(v=zs_total, col="green", lwd = 4,lty=2)
}

writeOut = function(file, mutD, tgene, zs_total, zs_full, optCorrRes = NA, tag = NA){
  #output
  outHeader = paste(c("tgene","totalSmp","mutSmp","optCorrSmp",
                      "totalReg","mutReg","optCorrReg",
                      "totalCorr","fullCorr:pval","permuCorr:pval","optCorr",
                      "selectGenes","selectSamples"),collapse="\t")
  
  optCorrSmp  <- optCorrReg <- outGene <- corr_opt <- corr_perm <- pvalF <- pvalP <- NA
  options(warn=-1)
  if (!is.na(optCorrRes)){
    resMut = optCorrRes$mutD
    corr_opt  = optCorrRes$corr_opt
    optCorrSmp = corr_opt$n
    corr_opt  = z2corr(corr_opt$zs)
    corr_perm = z2corr(optCorrRes$corr_perm$zs)
    pvalF = optCorrRes$pval_full
    pvalP = optCorrRes$pval_perm
    optCorrReg = length(optCorrRes$reg)
    resMut  = removeZeor(resMut)
    if (NROW(resMut) == 0){
      outSmpMap = c(paste("[",'', "]",sep=""), paste("[",'', "]",sep=""))
    
    }else{
    
    
      outGene = optCorrRes$reg
      outSmp  = optCorrRes$sample
      outSmpMap = vector(length=length(outGene))
      for (i in 1:nrow(resMut)){
        outSmpMap[i] = paste("[",paste(vapply(which(resMut[i,]!=0, arr.ind=T),
                                      function(x){outSmp[x]},'a'),collapse=";"),
                           "]",sep="")
      }
    }
  }else{
    resMut = removeZeor(mutD)
    outGene = rownames(mutD); outSmp = colnames(mutD)
    outSmpMap = vector(length=length(outGene))
    if(NROW(resMut) >1 & NCOL(resMut) > 1){
        for (i in 1:nrow(resMut)){
            outSmpMap[i] = paste("[", paste(vapply(which(resMut[i,]!=0, arr.ind=T),
                                function(x){outSmp[x]},'a'),collapse=";"),
                           "]",sep="")
        }
        
    }else{
      outSmpMap  = paste("[", paste(vapply(which(resMut !=0, arr.ind=T),
                                function(x){outSmp[x]},'a'),collapse=";"),
                           "]",sep="")
    }     
  }
  options(warn=0)
  
  if (! is.na(tag)) {
    tgene = paste(tag, tgene,sep=":")
  }
  
  outRecord = paste(c(tgene, zs_total$n, zs_full$n, optCorrSmp, 
                      nrow(mutD), nrow(removeZeor(mutD)), optCorrReg, 
                      format(z2corr(zs_total$zs),digits=4), 
                      paste(format(z2corr(zs_full$zs),digits=4), format.pval(pvalF),sep=":"),
                      paste(format(corr_perm,digits=4),format.pval(pvalP),sep=":"), 
                      format(corr_opt,digits=4),
                      paste(outGene, collapse=";"),
                      paste(outSmpMap,collapse=";") ), collapse="\t")
  
  write.table(outHeader, file=file, append=FALSE, col.names=F, row.names=F, sep="\t",quote=F)
  write.table(outRecord, file=file, append=TRUE, col.names=F, row.names=F, sep="\t",quote=F)
}
###func END-----


##load data----
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


##simple QC-------
if (is.null(removeZeor(mutD))){
  flag = FALSE
}else if (NROW(removeZeor(mutD)) < 2){
  flag = FALSE
}else if (NCOL(removeZeor(mutD)) < 10) {
  flag = FALSE
}else {
  flag = TRUE
}

##----
regExpD = subset(expD,select = names(tarExpD))[-1,]

if (NCOL(regExpD) == 1){ 
  sumExpD = regExpD 
}else{
  sumExpD = colSums(regExpD)
}
###------permute actual mutation started he

### all expression Mreg v.s Etar
zs_total = list(zs = fisherZ(cor(sumExpD, tarExpD)),
                  n = length(sumExpD) ) 

### full mutation matrix Mreg * Ereg v.s Etar
zs_full = calZscore(regExpD,mutD,tarExpD)

if (!flag) {
  #fail qc, exit 
  writeOut(output, mutD, tgene, zs_total, zs_full, tag="failQC")
  q(save = "no")
}

##select tolence------
### keeping the original method of selecting tol
nperm = numRandom
pvalCut = 0.1

if (typeTol == "flex") {
    tolVec = seq(-0.02,0.005,by=0.001)
    nTol = length(tolVec)
    
    iterTolSum = as.data.frame(matrix(NA, nrow=nTol,ncol=5))
    colnames(iterTolSum) = c("tol", "optcorr", "optSmpn",'optRegn', "pvalfull")
    
    for (i in 1:nTol ){
      tol = tolVec[i]
      resCorr_crt = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
      pval_full = resCorr_crt$zDiff_pval
      iterTolSum[i, ]  = c(tol = tol, optcorr=z2corr(resCorr_crt$corr_opt$zs),
                           optSmpn = resCorr_crt$corr_opt$n, 
                           optRegn = length(resCorr_crt$regs),
                           pvalfull = resCorr_crt$pval_full)
    }
        
    
    tol = iterTolSum[which(is.na(iterTolSum$pvalfull)),]
    
    # pvalue QC
    if (NROW(tol)  == 0 ){
        tol = iterTolSum[do.call(order, list(iterTolSum$pvalfull,iterTolSum$optRegn)),][1,1]
        resCorr_crt = corrOpt_binflip(mutD, regExpD, tarExpD, zs_full, tol = tol)
        writeOut(output, mutD, tgene,zs_total=zs_total, zs_full=zs_full, 
                 resCorr_crt, tag="failOptTol")
        q(save="no")
    }
    ##get tol
    tol = tol$tol[ which.max( vapply(tol$pvalfull, FUN=function(x){ifelse(is.na(x), 0, -log(x))},0.1) )[1] ]
    
}else if (typeTol == "fix") {
  tol = 0
}else{
  print(paste(ERR, "ttol < flex/fix> :", typeTol))
  q(save="no")
}

##random starting-----
print(paste("Finished Tol time :", format.difftime(difftime(Sys.time(), timeStart)), sep=""))

nIter = as.integer(numRandom)
pval_full_min <- pval_perm_min <- 1
    
if ( typeSelect == "max" ) {
  print("seletion using optimal")
  resRandInitOptCorr = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  

  for (i in 1:nIter){
    resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
    if(!is.na(resCorrOpt$pval_full) & 
           resCorrOpt$pval_full < pval_full_min){
        pval_full_min = resCorrOpt$pval_full
        resRandInitOptCorr = resCorrOpt
      }
  }
  resCorrOpt = resRandInitOptCorr
  resIterMutD = resRandInitOptCorr$mutD
}else{
  resIterMutD = matrix(0, nrow= nrow(mutD),ncol =ncol(mutD) )

  iterRansInitSum = as.data.frame(matrix(NA, nrow=nIter,ncol=5))
  colnames(iterRansInitSum) = c(paste("tol",tol,sep="_"), "optcorr",  
                                "optSmpn",'optRegn', "pvalfull")
  for (i in 1:nIter){
    resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, zs_full, tol = tol)  
    resIterMutD = resIterMutD + resCorrOpt$mutD
    iterRansInitSum[i, ]  = c(iter = paste(tol,i,sep="_"), 
                              optcorr = z2corr(resCorrOpt$corr_opt$zs),
                              optSmpn = resCorrOpt$corr_opt$n, 
                              optRegn = length(resCorrOpt$regs),
                              pvalfull = resCorrOpt$pval_full)
  }  
  if (typeSelect == 'all'){
    print("seletion uisng  all ")
    ncut = 0  
    resMut = apply(resIterMutD, c(1,2), function(x){ifelse(x >= ncut, 1, 0)})
    
    cntReg = resCorrOpt$corr_opt$n    
    final_corr = calZscore(expD=regExpD, mutD=resMut, tarExp=tarExpD)
    resCorrOpt <- list( mutD = resMut, corr_opt = final_corr,
                        regs = rownames(removeZeor(resMut)), sample = colnames(removeZeor(resMut)),
                        zs_full = zs_full, pval_full = zDiff(final_corr$zs, zs_full$zs, final_corr$n, zs_full$n))
    
  }else if (length(grep('[0-9.]+', typeSelect)) == 1 & length((grep('[a-z]+', typeSelect))) == 0){
    ncut = as.integer(typeSelect)
    print(paste("seletion using cutoff", ncut))
    
    resMut = apply(resIterMutD, c(1,2), function(x){ifelse(x >= ncut, 1, 0)})
    final_corr = calZscore(expD=regExpD, mutD=resMut, tarExp=tarExpD)
    resCorrOpt <- list( mutD = resMut, corr_opt = final_corr,
                        regs = rownames(removeZeor(resMut)), sample = colnames(removeZeor(resMut)),
                        zs_full = zs_full, pval_full = zDiff(final_corr$zs, zs_full$zs,final_corr$n,zs_full$n)) 
  }else{
    print( paste(ERR, "tinit should be <max> < all> or <0-100 numbers>", typeSelect))
#     q(save='no')
  }
}

##-------select and update resCorrOpt object

permuRes = t(sapply(seq(1,numRandom),FUN=function(x){
  mutD_perm  = permuMutD(mutD)
  print(t(which(mutD_perm>0, arr.ind=T)))
  return(optCorr(mutD_perm,regExpD,tarExpD, numRandom))
}))

save.image(file=paste(output,".Rdata",sep=""))
# 
# print(c(resCorrOpt$pval_full, unlist(permuRes[,3])))
# print(c(resCorrOpt$corr_opt$zs, unlist(vapply(permuRes[,1], FUN=function(x){x$zs},1))))
# print(c(resCorrOpt$corr_opt$n, unlist(vapply(permuRes[,1], FUN=function(x){x$n},1))))
# 
# 
# ##-------output
# output = paste(output , "_", tol, sep="")
# writeOut(output, mutD, tgene,zs_total=zs_total, zs_full=zs_full, resCorrOpt)
# write.table(file = paste(output,".fullMatrix", sep=""), x=removeZeor(resIterMutD), quote=F,col.names=T,sep="\t", row.names=T)
# 
# print(paste("[END]",tgene, " Running time :", format.difftime(difftime(Sys.time(), timeStart)), sep=""))
