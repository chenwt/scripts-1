install.packages("animation")
require(animation)
setwd("/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/May2014/doc/")
####RUN $crnsc/model/step3-4_greedyOptCorr.r for functions and get data!
library(animation)



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

corrDiff = function(z1,z2,n1,n2)  {
  if (n1 >=10 & n2 >=10){
    z = (z2 - z1) / sqrt( 1/(n2-3) + 1/(n1-3) ) 
    return(2*(1-pnorm(abs(z))))
  }else{
    return(NA)
  }
}

calCorr = function(expD, mutD, tarExp){
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

corrOpt_binflip = function(mutD, regExpD, tarExpD, corr_full = corr_full, tol = 0){
  ##get mutation indice
  allSmp = colnames(mutD)
  allReg = rownames(mutD)
  
  mutInd = which(mutD>0,arr.ind=T)
  numSmp = length(unique(mutInd[,2]))
  numReg = length(unique(mutInd[,1]))
  
  #QC
  if (numSmp < 10 | numReg < 2 ){
    return(list(corr_opt =NA, sample = numSmp , regs = numReg, 
                corrDiff_pval=NA, mutD=removeZeor(mutD)))
  }
  
  ##random initiate
  numMut = nrow(mutInd)
  mutInd = mutInd[sample(numMut),]
  mut_temp = mutD
  corr_prev = corr_full
  
  ##flipover forward selection
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
  
  ## pvalue with full
  corrDiffPval = corrDiff(corr_prev$zs, corr_full$zs,corr_prev$n,corr_full$n)
  smpOpt = names(colSums(mut_temp)[colSums(mut_temp)>0])
  numSmpOpt = length(smpOpt)
  regOpt = names(rowSums(mut_temp)[rowSums(mut_temp)>0])
  numRegOpt = length(regOpt)
  corr_opt = corr_prev
  return(list(corr_opt = corr_prev, corr_full=corr_full,
              pval_full = corrDiffPval,
              sample = smpOpt, regs =regOpt, 
              mutD=mut_temp))
  
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

doCorropt_perm = function(mutD, regExpD, tarExpD, corr_full, tol, nperm){
  resCorr_crt = corrOpt_binflip(mutD,regExpD,tarExpD,corr_full, tol = tol)
  corr_perm = doPermu(regExpD, resCorr_crt$mutD, tarExpD)
  pval_perm = max(length(corr_perm[abs(corr_perm) > abs(resCorr_crt$corr_opt$zs)])/nperm, 
                  1/nperm)
  resCorr_crt$corr_perm = list(zs=mean(corr_perm), n = resCorr_crt$corr_opt$n)
  resCorr_crt$pval_perm = pval_perm
  return(resCorr_crt)
}

plot_perm = function(corr_perm, resCorr, corr_full){
  hist(corr_perm, col="lightgray", xlim=c(-0.8,0.8),border="gray",
       xlab = "permut correlation", main = " correlation from permutation")
  abline(v=resCorr$corr_opt$zs,col="red",lwd=4)
  abline(v=corr_full$zs,col="orange",lwd=4)
  abline(v=corr_total, col="green", lwd = 4,lty=2)
}

##-----------

expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp"
mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp"

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
regExpD = subset(expD[-1,],select= names(tarExpD))
sumExpD = colSums(regExpD)
corr_total = cor(sumExpD,tarExpD)
sumExpD = colSums(tarExpD * mutD)

corr_full = calCorr(regExpD,mutD,tarExpD)
nperm = 1000
mutD_raw = mutD
mutD  = mutD_raw
mutInd = which(mutD>0,arr.ind=T)
numSmp = length(unique(mutInd[,2]))
numReg = length(unique(mutInd[,1]))
numMut = nrow(mutInd)
mutInd = mutInd[sample(numMut),]
mut_temp = mutD
corr_prev = corr_full
tol=0



require(gplots)

saveGIF({
  
  randfilpCorr = as.data.frame(matrix(NA, nrow=numMut, ncol = 3))
  randfilpCorr[1,] = c(1,z2corr(corr_prev$zs), "red" )
  for (i in 1:numMut){
    id_flip = mutInd[i,]
    mut_temp[id_flip[1], id_flip[2]] <- 0
    corr_temp = calCorr(regExpD,mut_temp,tarExpD)
    if (is.na(corr_temp$zs)) {break}
    if ( abs(corr_temp$zs) - abs(corr_prev$zs) < tol){
      mut_temp[id_flip[1], id_flip[2]] <- 1
      randfilpCorr[i+1, ] = c(i, z2corr(corr_temp$zs), "gray")  
      corr_temp = list(corr=0,n=0)
    }else{
      corr_prev = corr_temp
      randfilpCorr[i+1, ] = c(i, z2corr(corr_temp$zs), "red")  
      corr_temp = list(corr=0,n=0)
    }
    
    mutPlot = mut_temp[unique(mutInd[,1]),unique(mutInd[,2])]
    options(show.error.messages=F, dev='R_DEFAULT_DEVICE')
    par(font.axis = 2, mfrow=c(2,1))   
    heatmap.2(mutPlot, bg = "white",
              key=F,
              srtCol = 60, offsetRow = 0.01, offsetCol=0.01,
              col=colorRampPalette(c("white","red"))(50), 
              cexRow=1, cexCol=0.8, 
              sepwidth=c(0.01,0.01),
              colsep=0:ncol(mutPlot),
              rowsep=0:nrow(mutPlot),
              sepcolor="darkgray",
              trace="none", Rowv=F,dendrogram="none", Colv=F, 
              main=paste(tgene, "Flipping Over mutation", i,"\n", rownames(mutD)[id_flip[1]], ":", colnames(mutD)[id_flip[2]]))
    options(show.error.messages=T)
  }
}, interval = 0.2, ani.width = 1000,ani.height = 600)

saveGIF({
  par(font.axis = 2, font.main=2, cex.main=1.2, mar=c(3,10,4,2)) 
  for (i in 1: numMut){
    plot(x= 1:i, y = randfilpCorr[1:i,2], cex=1.2, cex.main= 2, font =2,
         col=ifelse(randfilpCorr[1:i,3] =="gray", "red","gray"), 
         xlab = "mutation", ylab = "correlation", main = paste("Greedy optimization of", tgene, "\ntol", tol),
         xlim=c(1,numMut),ylim = c(0,1),pch = 19)
  }
}, interval = 0.2, ani.width = 1000,ani.height = 600)


