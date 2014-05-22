CWD = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/"
setwd(CWD)

#test files------
expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_exp.temp"
mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1_regMut.temp"
figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/May2014/fig/"
output = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_1.tsv"
##test end------
usage = "Usage: Rscript step3-4_greedyOptCorr.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")

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

###---func----
require(gplots)
fisherZ = function(r) return(1/2*log((1 + r )/(1-r)))
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
  sumExp = colSums(expD * mutD)
  sumExp = sumExp[which(sumExp != 0)]
  numCol = length(sumExp)
  ##need further optimization  
  resCorr = fisherZ(cor(sumExp, tarExpD[names(sumExp)])) 
  return(list(zs=resCorr, n = numCol))
}

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

doPermu = function(regExpD, mutD, tarExpD, nperm = 1000){
  corr_perm = rep(0, nperm)
  for (i in 1:nperm){
    corr_perm[i]  = calCorr(expD=regExpD,mutD=permuMutD(mutD),tarExp=tarExpD)$zs
  }
  return(corr_perm)
}

plot_perm = function(corr_perm, resCorr){
  hist(corr_perm, col="lightgray", xlim=c(-0.8,0.8),border="gray",
       xlab = "permut correlation", main = " correlation from permutation")
  abline(v=resCorr$corr_opt$zs,col="red",lwd=4)
  abline(v=corr_full$zs,col="orange",lwd=4)
  abline(v=corr_total, col="green", lwd = 4,lty=2)
}

removeZeor = function(arr) return(arr[rowSums(arr)!=0, colSums(arr)!=0])

doCorropt_perm = function(mutD, regExpD, tarExpD, corr_full, tol, nperm){
  resCorr_crt = corrOpt_binflip(mutD,regExpD,tarExpD,corr_full, tol = tol)
  corr_perm = doPermu(regExpD, resCorr_crt$mutD, tarExpD)
  pval_perm = max(length(corr_perm[abs(corr_perm) > abs(resCorr_crt$corr_opt$zs)])/nperm, 
                  1/nperm)
  resCorr_crt$corr_perm = list(zs=mean(corr_perm), n = resCorr_crt$corr_opt$n)
  resCorr_crt$pval_perm = pval_perm
  return(resCorr_crt)
}


###----loadData-----

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
tol =  -0.001

##--tselmax plot---
nIter = 100
pval_full_min <- pval_perm_min <- 1
for (i in 1:nIter){
  resCorrOpt = doCorropt_perm(mutD, regExpD, tarExpD, corr_full, tol, nperm)
  if(!is.na(resCorrOpt$pval_full) & !is.na(resCorrOpt$pval_perm) & 
       resCorrOpt$pval_full < pval_full_min & resCorrOpt$pval_perm < pval_perm_min){
    pval_full_min = resCorrOpt$pval_full
    pval_perm_min = resCorrOpt$pval_perm
    resRandInitOptCorr = resCorrOpt
  }
}

###---compute selected mutation 
smpMutAll = colnames(removeZeor(mutD))
smpMutSelL = resRandInitOptCorr$sample
smpMutNonSelL = setdiff(, smpMutSelL)
smpNonMutL = setdiff(names(sumExpD), colnames(removeZeor(mutD)))
regSel = resRandInitOptCorr$regs
mutDSel = resRandInitOptCorr$mutD
sumExp = colSums(regExpD * mutDSel)
sumExp = tarExpD
drCerRNA_selMut_nonMut <- 
  sapply(regSel,function(x){
    smpMutSel <- colnames(mutDSel[x,] != 0)
    t.test(sumExp[smpMutSelL], sumExp[smpNonMutL])
    })

drCerRNA_selMut_nonMut

tgene_selMut_nonMut <- t.test(tarExpD[smpNonMutL], tarExpD[smpSelMutL])

diff_selMut_nonMut$p.value

outMutD = removeZeor(resRandInitOptCorr$mutD)
plotMutD3D = function(outMutD){
  require("plot3D")
  outMutD = outMutD[order(rowSums(outMutD),decreasing=T),]
  z = outMutD 
  z = z[,order(colSums(z),decreasing=T)]
  zz = apply(z, c(1,2),function(x){ifelse(x!=0,1,0)})
  z = z[,do.call(order, lapply(1:nrow(z), function(i) -z[i, ]))]
  myColVar = t(apply(z, 1, function(x) {x  * as.numeric(as.factor(colnames(z)))}))
  hist3D( x = seq(0, 1, length.out = nrow(z)), zlim = c(0,1.5),
          y = seq(0, 1, length.out = ncol(z)), z, 
          colvar =myColVar  , phi = 40, theta = 40,
          xlab = "ceRNA drivers", ylab ="selected samples", 
          col = NULL, NAcol = "white", border = NA, facets = TRUE,
          colkey = NULL, image = F, contour = F,
          panel.first = NULL, clim = NULL, clab = NULL, bty = "b",
          lighting = T, shade = F, ltheta = -135, lphi = 0,
          space = 0.2, add = F, plot = TRUE)
  # text3d( x = seq(0, 1, length.out = nrow(z)), zlim = c(0,1.5),
  #         y = seq(0, 1, length.out = ncol(z)),  text=fileName,adj = 0.1, 
  #        color=setRamp(6), family="serif", font=5, cex=1)
  # heatmap.2(z,trace="none",col=c("white",'red'))
}

plotMutD3D(t(outMutD))


##---tselAll plot---
pdf(paste(figd, "/plot_step3-4_greedyOptCorrResult_",tgene, "_" ,CDT, "_heatmap.pdf",sep=""))
par(mar=c(10,6,4,6),mgp=c(10,0.2,0))
# heatmap.2(mutD, col=colorRampPalette(c("white","red"))(nIter), trace="none",Colv=F)
# heatmap.2(resIterMutD, col=colorRampPalette(c("white","red"))(nIter), trace="none",Colv=F)
heatmap.2(removeZeor(resIterMutD), col=colorRampPalette(c("white","red"))(nIter), trace="none", Rowv=F,dendrogram="none", Colv=F, main=paste("Random init\n target:", tgene, "\nIter(n):",nIter))
ncut = 30
resMut = apply(resIterMutD, c(1,2), function(x){ifelse(x > ncut, 1, 0)})
heatmap.2(removeZeor(resMut), col=colorRampPalette(c("white","red"))(nIter), trace="none",Rowv=F,dendrogram="none",Colv=F,main=paste("final result from Random init\n target:", tgene, "\ncutoff:",ncut,"/",nIter))
dev.off()

dim(removeZeor(resCorrOpt$mutD))
dim(removeZeor(resIterMutD))
dim(removeZeor(resMut))

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
write.table(outRecord, file=output, append=T, col.names=F, row.names=F, sep="\t",quote=F)
# 
write.table(removeZeor(resMut), file=paste(output,".mutmatrix",sep=""), quote=F, sep="\t")
# 
par(mfrow=c(1,1),mar=c(6,2,2,4))
# plot(hclust(dist(t(outD))))
heatmap.2(removeZeor(resMut), col=colorRampPalette(c("white","red"))(nIter), trace="none")
dev.off()

# install.packages("biclust")
# require("biclust")
# drawHeatmap(outD, biclust(outD,method="BCBimax"),2)

# outD = outD[rowSums(outD) >10, colSums(outD) > 10]
# plot_perm(vapply(corr_perm,z2corr,0.2), resCorrOpt)

## calculate final output
# outTarExpD  = tarExpD[resCorrOpt$sample ]
# outTarExpD  = outTarExpD[order(outTarExpD)]
# outMutD = resCorrOpt$mutD
# outMutD = outD[rowSums(outMutD)!=0, colSums(outMutD) != 0 ]
# outD = rbind(CEP55=outTarExpD,subset(regExpD,select=names(outTarExpD),subset=rownames(regExpD) %in% names(resCorrOpt$cerna)))
# outMutD = resCorrOpt$mutD
# outMutD = outMutD[rowSums(outMutD) !=0,colSums(outMutD) !=0]
# heatmap.2(outD,col=bluered(250),trace="none",scale="row",Rowv=F,dendrogram="none")
# heatmap.2(expD,col=bluered(250),trace="none",scale="row",Rowv=F,dendrogram="none")
# dev.off()
# save.image("step3-4_gdOptCorr.rda")

install.packages("plot3D")
plotMutD3D = function(outMutD){
  require("plot3D")
  outMutD = outMutD[order(rowSums(outMutD),decreasing=T),]
  z = outMutD 
  z = z[,order(colSums(z),decreasing=T)]
  zz = apply(z, c(1,2),function(x){ifelse(x!=0,1,0)})
  z = z[,do.call(order, lapply(1:nrow(z), function(i) -z[i, ]))]
  myColVar = t(apply(z, 1, function(x) {x  * as.numeric(as.factor(colnames(z)))}))
  hist3D( x = seq(0, 1, length.out = nrow(z)), zlim = c(0,1.5),
          y = seq(0, 1, length.out = ncol(z)), z, 
          colvar =myColVar  , phi = 40, theta = 40,
          xlab = "ceRNA driver", ylab ="selected samples", 
          col = NULL, NAcol = "white", border = NA, facets = TRUE,
          colkey = NULL, image = F, contour = F,
          panel.first = NULL, clim = NULL, clab = NULL, bty = "b",
          lighting = T, shade = F, ltheta = -135, lphi = 0,ticktype = "detailed",
          space = 0.2, add = F, plot = TRUE)
  # text3d( x = seq(0, 1, length.out = nrow(z)), zlim = c(0,1.5),
  #         y = seq(0, 1, length.out = ncol(z)),  text=fileName,adj = 0.1, 
  #        color=setRamp(6), family="serif", font=5, cex=1)
  # heatmap.2(z,trace="none",col=c("white",'red'))
}
plotMutD3D(removeZeor(resIterMutD)/50)

pdf(paste(figd, "/plot_step3-4_greedyOptCorr_",CDT, "_3D.pdf",sep=""))
par(mfrow=c(1,2))
plotMutD3D(removeZeor(mutD))

plotMutD3D(removeZeor(resMut))

mutSumD = removeZeor(mutD + resMut)
mutSumBinD = apply(mutSumD,c(1,2),function(x){ifelse(x!=0,1,0)
                                              
                                              plotMutD3D(mutSumD)
                                              par(mfrow=c(1,1))
                                              
                                              dev.off()
                                              
                                              
                                              
                                              ###=============
                                              
                                              library(rgl)
                                              
                                              
                                              hist3d<-function(x,y=NULL,nclass="auto",alpha=1,col="#ff0000",scale=10)
                                              {
                                                save <- par3d(skipRedraw=TRUE)
                                                on.exit(par3d(save))
                                                xy <- xy.coords(x,y)
                                                x <- xy$x
                                                y <- xy$y
                                                n<-length(x)
                                                if (nclass == "auto") { nclass<-ceiling(sqrt(nclass.Sturges(x))) }
                                                breaks.x <- seq(min(x),max(x),length=(nclass+1))
                                                breaks.y <- seq(min(y),max(y),length=(nclass+1))
                                                z<-matrix(0,(nclass),(nclass))
                                                for (i in 1:nclass) 
                                                {
                                                  for (j in 1:nclass) 
                                                  {
                                                    z[i, j] <- (1/n)*sum(x < breaks.x[i+1] & y < breaks.y[j+1] & 
                                                                           x >= breaks.x[i] & y >= breaks.y[j])
                                                    binplot.3d(c(breaks.x[i],breaks.x[i+1]),c(breaks.y[j],breaks.y[j+1]),
                                                               scale*z[i,j],alpha=alpha,topcol=col)
                                                  }
                                                }
                                              }
                                              
                                              rgl.bg(col="#cccccc")
                                              x<-rnorm(250)
                                              y<-rnorm(250)
                                              hist3d(x, y, alpha=0.8, nclass=10, scale=30)
                                              axes3d(c('x','y','z'))
                                              title3d('','','xlab','ylab','zlab')plot
                                              
                                              
                                              library(lattice) 
                                              library(latticeExtra) 
                                              ?panel.3dbars
                                              