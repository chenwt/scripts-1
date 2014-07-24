#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
###J.HE
###----
# get driver mutations for TF regulator of target genes
# input : 1 expression 
#         2 mutation 
# output :
###----

rm(list=ls())
usage = "Usage: Rscript step3-4_greedyOpt.r  -exp <expression file from python> -mut <mutation file from python>"
ERR = "ERROR:"
CDT = paste(unlist(strsplit(system('date',intern=T)," "))[c(2,4,7)],collapse="-")
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


#----getting command line parameters
args = getArgs()
if(length(args) < 1 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop(paste(error,"wrong input parameter!"))
}

setwd(system("pwd",intern=T))
dataFile    = args[['data']]
nrandStart  = as.integer(args['nrandStart'])
nrandk      = as.integer(args['nrandk'])


print("###------------------------")
print(paste("data file: ", class(dataFile), dataFile))
print(paste("nrandStart: ", nrandStart))
print(paste("nrandk: ", nrandk))
print(paste("output: ", output))

#---init

# ###test files------
# nrandk = 100
# nrandStart = 100
# dataFile = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/test/output.tfmut.test_PTPLB.temp"
# output = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/result/07152014/tfMut/test/selMS_PTPLB.txt"
# figd = "/Volumes//ifs/data/c2b2/ac_lab/jh3283/projFocus/report/Jul2014/fig/"
# 
# # ##test end------

###func-----
getMutSmp = function(mutD){
  mutSmps = which(apply(mutD,2,function(x){any(x>0) == T}) == T, arr.ind=T)
  return(names(mutSmps))
}

calKS = function(x,y) { options(warn=-1); ts_act = ks.test(x=x, y=y, alternative="two.sided",exact= NULL); options(warn=1); return(ts_act)}

optMutSample = function(mutD, expD){
  ######
  # one round of greedy with one random start 
  # output the optimized sample and statistics
  #####
  isTsBetter = function(ts_prev, ts_crt ){return (ts_prev$statistic <= ts_crt$statistic) }
  
  ##get mutation indice
  allSmp = colnames(mutD)
  allReg = rownames(mutD)
  
  mutInd = which(mutD>0,arr.ind=T)
  numSmp = length(unique(mutInd[,2]))
  numReg = length(unique(mutInd[,1]))
  
  
  ms_act = getMutSmp(mutD)
  noms_act = setdiff(allSmp, ms_act)
  
  ts_act = calKS(x=expD[ms_act], y=expD[noms_act])
  
  
  ##random initiate
  randInd = sample(length(ms_act))
  ms_prev = ms_act; ts_prev = ts_act
  for (i in randInd){
    ms_crt = ms_prev[-i]
    noms_crt = setdiff(allSmp, ms_crt)
    ts_crt = calKS(x=expD[ms_crt], y= expD[noms_crt])
    if (isTsBetter(ts_prev, ts_crt)){
      ms_prev = ms_crt;  ## use new samples, else keep. 
      ts_prev = ts_crt; ts_prev = ts_crt
    }
  }
  ts_opt = ts_prev; ms_opt = ms_prev
  return(list(ms=ms_opt, ts = ts_opt))
}

permuLables = function(mutD){
  ### permute row and column labels for input matrix
  colname_orig = colnames(mutD)
  rowname_orig = row.names(mutD)
  mutD_perm = mutD
  mutD_perm = mutD_perm[sample(nrow(mutD_perm)),]
  mutD_perm = mutD_perm[,sample(ncol(mutD_perm))]
  colnames(mutD_perm) <- colname_orig
  rownames(mutD_perm) <- rowname_orig 
  return(mutD_perm)
}

rand_SelectKMut = function(mutD, k) {
  ######
  # random select k mutations from original n mutation matrix
  ######
  allSmp = colnames(mutD);  allReg = rownames(mutD)
  cnts = ncol(mutD); cntg = nrow(mutD)
  
  mutInd = which(mutD>0,arr.ind=T)  
  mutInd = mutInd[sample(1:nrow(mutInd),k),]
  
  mutD[1:cntg, 1:cnts] = 0
  for (i in 1:k){
    mutD[mutInd[i,1], mutInd[i,2]] <- 1
  } 
  
  return(mutD)
}


randK = function (optMSobj, mutD){
  #########
  # get test statistics for randome choose k mutation from original mutation matrix
  # k is the same number of optimized mutations
  #########
  ms_randk  = getMutSmp(rand_SelectKMut(mutD=mutD, k=NROW(which(mutD[optMSobj$ms] >0,arr.ind=T))) )
  allsample = colnames(mutD)
  ts_randk = calKS(expD[ms_randk],expD[setdiff(allsample, ms_randk)])
  return(ts_randk)
}

plot_perm = function(corr_perm, resCorr, corr_full){
  hist(corr_perm, col="lightgray", xlim=c(-0.8,0.8),border="gray",
       xlab = "permut correlation", main = " correlation from permutation")
  abline(v=resCorr$corr_opt$zs,col="red",lwd=4)
  abline(v=corr_full$zs,col="orange",lwd=4)
  abline(v=corr_total, col="green", lwd = 4,lty=2)
}

getGeneSmps = function(mutD, mutSample ){
  allgene = rownames(mutD); allsample = colnames(mutD)
  mutIndex = which(mutD[,mutSample] >0, arr.ind=T)
  cnt = NROW(mutIndex)
  muts = cbind(allgene[mutIndex[,1]], allsample[mutIndex[,2]])
  colnames(muts) <- c("gene", "sample")
  res = NA
  for (g in unique(muts[,1])){
    gsample = paste("[", g, ":", paste(muts[which(muts[,1] == g),2], collapse=","), "]",sep="")
    res = c(res, gsample)
  }
  return(paste( res[-1], collapse= ";"))
}

Stouffer.test <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p/2) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(c(Z = Z, p.value = p.val))
}

###func END-----


###load data----
dataDF = read.table(dataFile, sep="\t",header = T, stringsAsFactors=F); rownames(dataDF) = dataDF[,1] ; dataDF = data.frame(as.matrix(dataDF[,-1]))

print(dim(dataDF))
expD = dataDF[1,]; colnames(expD) = colnames(dataDF); expD = unlist(expD)
mutD = dataDF[-1,]
tgene = rownames(dataDF)[1]
regs = rownames(mutD)


###main----

## init--
res_optMS_prev = optMutSample(mutD,expD)

## run multiple times --
for (i in 1:nrandStart){
  print(paste("Random starting ", i))
  res_optMS_crt = optMutSample(mutD,expD) ## everytime this run, 'sample' will run random starting 
  
  if(res_optMS_crt$ts$statistic > res_optMS_prev$ts$statistic){
    res_optMS_prev = res_optMS_crt
    print(paste(length(res_optMS_prev$ms), paste(unlist(res_optMS_prev$ts)[1:2],collapse=" "), sep=" "))
  }
}

## run choose k from all mutations
print("## run choose k from all mutations ")
res_optMS_crt <- res_optMS_prev; rm(res_optMS_prev)
res_optMS_randK = rep(NA, 2); for ( irand in 1:nrandk) res_optMS_randK = rbind(res_optMS_randK ,unlist(randK(res_optMS_crt, mutD) )[1:2])  
pval_randK = max(length(which(as.numeric(res_optMS_crt$ts$statistic) < as.numeric(res_optMS_randK[,1][-1]), arr.ind=T)) / nrandk, 1/nrandk)
res_optMS_crt$pval_randK = pval_randK; rm(pval_randK)

###output----
##output selected sample mutations vector 
tarRegMut = data.frame(matrix(0, ncol= NCOL(mutD))) ; colnames(tarRegMut) <- colnames(mutD); 
tarRegMut[res_optMS_crt$ms] <- unlist(colSums(mutD[,res_optMS_crt$ms]))
row.names(tarRegMut) <- tgene

output1 = paste(output, tgene, nrandStart, nrandk, "mutSampleVector", sep="_")
write.table(x=tarRegMut, file=output1, sep="\t", quote = F, row.names=T,col.names=T )


##output summaries
ms_act = getMutSmp(mutD); noms_act = setdiff(colnames(mutD), ms_act); ts_act = calKS(x=expD[ms_act], y=expD[noms_act])
cntRegs = length(which(rowSums(mutD) > 0 , arr.ind=T))
cntSmps = length(which(colSums(mutD) > 0 , arr.ind=T))
out_original = c(cntRegs, cntSmps, unlist(ts_act)[1:2])
  
cntOptRegs = length( which(rowSums(mutD[,res_optMS_crt$ms]) > 0 , arr.ind=T))
out_opt = c(cntOptRegs, length(res_optMS_crt$ms), unlist(res_optMS_crt$ts)[1:2], res_optMS_crt$pval_randK)

p <- c(ts_act$p.value, res_optMS_crt$ts$p.value)
w <- c(cntSmps, length(res_optMS_crt$ms))
pval_acVSopt = Stouffer.test(p = p, w = w) 
out_compare = pval_acVSopt

out_param = c(nrandStart, nrandk)

out_gSmps = getGeneSmps(mutD, res_optMS_crt$ms)

out_header = c("tgene", "mutRegs_act", "mutSmps_act", "D_act", "pval_act",
               "mutRegs_opt", "mutSmps_opt", "D_opt", "pval_opt", "pval_randK",
               "Z_actVopt", "pval_actVopt",
               "regSmps_opt",
               "numRandStart", "numRandK")
out = c(tgene, out_original, out_opt, out_compare,out_gSmps, out_param)

output = paste( output, tgene, nrandStart, nrandk, sep="_")
write.table(x=rbind(out_header, out), file=output, sep="\t", quote = F, row.names=F,col.names=F )

print(paste("[END] ",tgene, " Running time :", format.difftime(difftime(Sys.time(), timeStart)), sep=""))
