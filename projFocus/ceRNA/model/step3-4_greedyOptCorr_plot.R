# CWD = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/"
# setwd(CWD)

##test files------
# expfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_BRD3_exp.temp"
# mutfile = "/Users/jh3283/HOME/DATA/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_BRD3_regMut.temp"
# figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/May2014/fig/"
# output = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/05012014/sigMut/test/test_optCorr_05062014_BRD3.tsv"
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
# #----getting command line parameters
args = getArgs()
if(length(args) < 1 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop(paste(error,"wrong input parameter!"))
}

setwd(system("pwd",intern=T))
expfile     = args['exp']
mutfile     = args['mut']
output      = args['output']
figd        = paste(rootd, "/DATA/projFocus/report/May2014/fig/", sep="")
print(paste("inputfile",expfile, mutfile))
print(paste("outputfile",output))
#---init

###---func
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
    
#     plot(x=1:numMut, ylim = c(0,1),type="n",xlab="number of Mutation",ylab="Correlation",main = paste("Greedy optimization for correlation\ntol =",tol))
#     points(1, z2corr(corr_prev$zs),col="red", pch = 16)
#     
    for (i in 1:numMut){
      id_flip = mutInd[i,]
      mut_temp[id_flip[1], id_flip[2]] <- 0
      corr_temp = calCorr(regExpD,mut_temp,tarExpD)
      if (is.na(corr_temp$zs)) {break}
      if ( abs(corr_temp$zs) - abs(corr_prev$zs) < tol){
        mut_temp[id_flip[1], id_flip[2]] <- 1
        points(i+1, z2corr(corr_temp$zs), col="gray",pch=16)         
        corr_temp = list(corr=0,n=0)
      }else{
        points(i+1, z2corr(corr_temp$zs), col="red",pch=16)        
        corr_prev = corr_temp
        corr_temp = list(corr=0,n=0)
      }
      heatmap.2(mut_temp, 
                col=colorRampPalette(c("white","red"))(50), 
                trace="none", Rowv=F,dendrogram="none", Colv=F, main=paste("Random init\n target:", tgene, "\nIter(n):",nIter))
    }
    corrDiffPval = corrDiff(corr_prev$zs, corr_full$zs,corr_prev$n,corr_full$n)
    
#     text(x=3,y=0.85,labels=paste("opt vs.full, p-val:", format(corrDiffPval,digits=4)), font = 2, pos=4,cex=0.8)    
    smpOpt = colSums(mut_temp)[colSums(mut_temp)>0]
    numSmpOpt = length(smpOpt)
    regOpt = rowSums(mut_temp)[rowSums(mut_temp)>0]
    numRegOpt = length(regOpt)
    corr_opt = corr_prev
    corr_opt$corr <- (exp(2 * corr_prev$zs) -1)/(exp(2 * corr_prev$zs) +1)
#     text(x = 3, y = 0.95, labels=paste("Opt Smp:",numSmpOpt,"/",numSmp),font=2,pos=4,cex=0.8)
#     text(x = 3, y = 0.9, labels=paste("Opt Gen:",numRegOpt,"/", numReg),font=2,pos=4,cex=0.8)
    return(list(corr_opt =corr_opt, sample = smpOpt, cerna=regOpt,corrDiff_pval=corrDiffPval, mutD=mut_temp))
  }
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
regExpD = subset(expD[-1,],select= names(tarExpD))
sumExpD = colSums(regExpD)
corr_total = cor(sumExpD,tarExpD)
sumExpD = colSums(tarExpD * mutD)


corr_full = calCorr(regExpD,mutD,tarExpD)
nperm = 1000

pdf(paste(figd, "/plot_step3-4_greedyOptCorr_",CDT, "_tol.pdf",sep=""),width=10,height=10)

par(mfrow=c(4,3),mar=c(2,2,4,0))
tolVec = seq(-0.01, 0.005,by=0.001)
nTol = length(tolVec)
iterTolSum = as.data.frame(matrix(NA, nrow=nTol,ncol=6))
colnames(iterTolSum) = c("tol", "optcorr", "optSmpn","optRegn","pvalPerm", "pvalfull")
for (i in 1:nTol ){
  tol = tolVec[i]
  resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD,corr_full, tol = tol)
  corr_perm = doPermu(regExpD, resCorrOpt$mutD, tarExpD)
  pval_perm = max(length(corr_perm[corr_perm > resCorrOpt$corr_opt$zs])/nperm, 0.00001)
  text(x=3,y=0.8,labels=paste("opt vs. perm p-val:", format(pval_perm,digits=4)), font = 2, pos=4,cex=0.8)
  pval_full = resCorrOpt$corrDiff_pval
  iterTolSum[i, ]  = c(tol = tol, optcorr=z2corr(resCorrOpt$corr_opt$zs), optSmpn = round(resCorrOpt$corr_opt$n),optRegn = length(resCorrOpt$cerna), pvalPerm = pval_perm, pvalfull=pval_full)
}

dev.off()
out = paste(output, ".iterTols.summary", sep="")
write.table(iterTolSum, out, col.names=F,sep="\t",quote=F)

tol = iterTolSum[which(iterTolSum[,4] < 0.001 & iterTolSum[,5] < 0.001),]
if (length(tol[,1]) >0) {
  tol = tol[which.max(tol[,2]),1]
}else{
  print("greedy didn't coverge!!")
  q(save="no")
}

nIter = 100
par(mfrow=c(3,4),mar=c(1,1,2,0))
resIterMutD = matrix(0, nrow=nrow(resCorrOpt$mutD),ncol = ncol(resCorrOpt$mutD))
iterRansInitSum = as.data.frame(matrix(NA, nrow=nTol,ncol=5))
colnames(iterRansInitSum) = c("#inter", "optcorr", "optSmpn","pvalPerm", "pvalfull")
pvalCut = 0.001; corr_max = 0.0
for (i in 1:nIter ){
  resCorrOpt = corrOpt_binflip(mutD,regExpD,tarExpD,corr_full, tol = tol)
  corr_perm = doPermu(regExpD, resCorrOpt$mutD, tarExpD)
  pval_perm = max(length(corr_perm[corr_perm > resCorrOpt$corr_opt$zs])/nperm, 0.00001)
  pval_full = resCorrOpt$corrDiff_pval
  if(!is.na(pval_full) & !is.na(pval_perm) & pval_full > pvalCut & pval_perm > pvalCut & z2corr(resCorrOpt$corr_opt$zs) > corr_max){
    corr_max =  z2corr(resCorrOpt$corr_opt$zs)
    print(paste(i,"corr_max", corr_max))
    resMutD = resCorrOpt$muD
  }
    resIterMutD = resIterMutD + resCorrOpt$mutD
    iterRansInitSum[i, ]  = c(i, optcorr=z2corr(resCorrOpt$corr_opt$zs), optSmpn = round(resCorrOpt$corr_opt$n), pvalPerm = pval_perm, pvalfull=pval_full)

}


out = paste(output, ".randomInti.summary", sep="")
write.table(iterRansInitSum, out, row.names=F,sep="\t",quote=F)

pdf(paste(figd, "/plot_step3-4_greedyOptCorr_",tgene, "_" ,CDT, "_heatmap.pdf",sep=""))
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
