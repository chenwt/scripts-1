CWD = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/"
setwd(CWD)
# keyRegSumfile="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/run-Apr-1-2014/summary/target_ceRNADriver_0.01"
# expfile="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/exp/brca_exp_l3_tumor_Mar-21-2014.matrix_Mar-26-2014.voomNormed.matrix"
# gslistfile="/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/gslist/gslist_CnvMethSomFree.10smapMore.deg_20140430.txt.hasReg.list"
# mutfile="/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/tcgal2som/brca_somlevel2_byGene.matrix.inExpSmp.matrix.nonzero"

expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/CEP55_exp"
mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/CEP55_regMut"
###---func

fisherZ = function(r) return(1/2*log((1 + r )/(1-r)))
corrDiff = function(r1,r2,n1,n2)  {
  z = (r2 - r1) / sqrt( 1/(n2-3) + 1/(n1-3) ) 
  return(2*(1-pnorm(abs(z))))
}
calCorr = function(expD, mutD, tarExp){
  sumExp = colSums(expD * mutD)
  sumExp = sumExp[which(sumExp != 0)]
  numCol = length(sumExp)
  resCorr = fisherZ(cor(sumExp, tarExpD[names(sumExp)])) 
  return(list(corr=resCorr, n = numCol))
}


###----func



require(useful)
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
cor(sumExpD,tarExpD)


corrOpt_binflip = function(mutD, regExpD, tarExpD, tol = -0.005){
  mutInd = which(mutD>0,arr.ind=T)
  numMut = nrow(mutInd)
  numSmp = length(unique(mutInd[,2]))
  numReg = length(unique(mutInd[,1]))
  mut_temp = mutD
  corr_full = calCorr(regExpD,mutD,tarExpD)
  corr_prev = corr_full
  print(c("full","full",corr_prev$corr, "", ""))
  options(digits=3)
  plot(x=1:numMut, ylim = c(0,1),type="n",xlab="iteration",ylab="Correlation",main = paste("Greedy optimization for correlation\ntol =",tol))
  points(1, corr_prev$corr,col="red", pch = 16)
  for (i in 1:numMut){
    id_flip = mutInd[i,]
    mut_temp[id_flip[1], id_flip[2]] <- 0
    corr_temp = calCorr(regExpD,mut_temp,tarExpD)
    print(c(rownames(mutD)[id_flip[1]],colnames(mutD)[id_flip[2]],corr_prev$corr, corr_temp$corr))
    if ( abs(corr_temp$corr) - abs(corr_prev$corr) < tol){
      mut_temp[id_flip[1], id_flip[2]] <- 1
      points(i+1, corr_temp$corr, col="gray",pch=16) 
      corr_temp = list(corr=0,n=0)
     }else{
      points(i+1, corr_temp$corr, col="red",pch=16)
      corr_prev = corr_temp
      corr_temp = list(corr=0,n=0)
    }
  }
  corrDiff_opt = corrDiff(corr_prev$corr, corr_full$corr,corr_prev$n,corr_full$n)
  smpOpt = colSums(mut_temp)[colSums(mut_temp)>0]
  numSmpOpt = length(smpOpt)
  regOpt = rowSums(mut_temp)[rowSums(mut_temp)>0]
  numRegOpt = length(regOpt)
  corr_Opt = corr_prev
  text(x = 3, y = 0.95, labels=paste("Opt Smp:",numSmpOpt,"/",numSmp),font=2,pos=4)
  text(x = 3, y = 0.9, labels=paste("Opt Gen:",numRegOpt,"/", numReg),font=2,pos=4)
  return(list(corr_opt =corr_prev, sample = smpOpt, cerna=regOpt,corr_full= corr_full,mutD=mut_temp))
}

resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, tol = 0.001)
resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, tol = 0)
resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, tol = -0.001)
resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, tol = -0.002)
resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, tol = -0.003)

resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD, tol = -0.004)

doPermu = function(regExpD, mutD, tarExpD, nperm = 1000){
  permuMutD = function(mutD){
    colname_orig = colnames(mutD)
    rowname_orig = row.names(mutD)
    mutD_perm = mutD
    mutD_perm = mutD_perm[sample(nrow(mutD_perm)),]
    mutD_perm = mutD_perm[,sample(ncol(mutD_perm))]
    colnames(mutD_perm) <- colname_orig
    rownames(mutD_perm) <- rowname_orig 
    return(mutD_perm)
  }
  corr_perm = rep(0, nperm)
  for (i in 1:nperm){
    corr_perm[i]  = calCorr(expD=regExpD,mutD=permuMutD(mutD),tarExp=tarExpD)$corr
  }
  return(corr_perm)
}



nperm = 1000
corr_perm = doPermu(regExpD, resCorrOpt$mutD, tarExpD)

plot_perm = function(corr_perm, resCorr){
  hist(corr_perm, col="lightgray", xlim=c(-0.8,0.8),border="gray",
       xlab = "permut correlation", main = " correlation from permutation")
  abline(v=resCorrOpt$corr_opt$corr,col="red",lwd=4)
  abline(v=resCorrOpt$corr_full$corr,col="orange",lwd=4)
  abline(v=corr_total, col="green", lwd = 4,lty=2)
}
corrDiff_perm = corrDiff(resCorrOpt$corr_opt$corr, resCorrOpt$corr_full$corr,resCorrOpt$corr_opt$n, resCorrOpt$corr_full$n)
colSums(mut_temp)[colSums(mut_temp)>0]
rowSums(mut_temp)[rowSums(mut_temp)>0]


outTarExpD  = tarExpD[resCorrOpt$sample ]
outTarExpD  = outTarExpD[order(outTarExpD)]

outD = rbind(CEP55=outTarExpD,subset(regExpD,select=names(outTarExpD),subset=rownames(regExpD) %in% names(resCorrOpt$cerna)))
heatmap.2(outD,col=bluered(250),trace="none",scale="row",Rowv=F,dendrogram="none")
heatmap.2(expD,col=bluered(250),trace="none",scale="row",Rowv=F,dendrogram="none")


# require(gplots)
# heatmap.2(regExpD, col=bluered(100),trace="none",scale='row',Colv=F)
# heatmap.2(rbind(blank=rep(0,length(tarExpD)),target=tarExpD), col=bluered(100),trace="none",scale='row', Rowv= F, colv = F)
