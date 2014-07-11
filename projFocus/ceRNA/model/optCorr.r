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

permuRes = NA
sapply(seq(1,numRandom),FUN=function(x){
  mutD_perm  = permuMutD(mutD)
  print(t(which(mutD_perm>0, arr.ind=T)))
  permuRes = c(permuRes, optCorr(mutD_perm,regExpD,tarExpD, numRandom))
})

print(c(as.vector(resCorrOpt$pval_full), vapply(permuRes,FUN=function(x){x$pval_full}, 0.11)))

