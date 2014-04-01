#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: < expression matrix (first row as the target, others as regulator) > 


example = "Example: /ifs/home/c2b2/ac_lab/jh3283/tools/R/R_current/bin/Rscript 
    /ifs/home/c2b2/ac_lab/jh3283/scripts/projFocus/ceRNA/model/step2-2_regKeyRegulators.r
   /ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/test/PTENreg_exp_tumor"

usage = "Usage: Rscript step2-2_regKeyRegulators.r  <exp.matrix>"
ERR = "ERROR here: "
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

#---funcs
regfit      = function(X, y, group){
   if (setequal(colnames(X), names(group))){
    rss   = NULL
    a     = seq(0.01,1,0.1)
    for (alpha in a){
      grpfit            = grpreg(X,y,group=group,penalty="grLasso",
                                 alpha = alpha,eps=0.05,max.iter=10000)
      grpfit_selected   = select(grpfit,criterion="AIC")
      rss               = c(rss,sum(y - grpfit_selected$beta[1] -
                                      X %*% as.matrix(grpfit_selected$beta[-1])))
    }
    alpha = a[which.min(rss)[1]]
    #fit model
    grpfit              = grpreg(X, y, group=group,penalty="grLasso",
                                 alpha = alpha,eps=0.005,max.iter=5000)
    grpfit_selected     = select(grpfit,criterion="AIC")   
    beta                = grpfit_selected$beta
    residual            = y - beta[1] - X %*% as.matrix(beta[-1])
    residual_sum_squre  = sum((residual)^2)
    reg_sum_square      = sum((X %*% as.matrix(beta[-1]) + beta[1] - mean(y))^2)
    total_sum_square    = reg_sum_square + residual_sum_squre
    r_square            = reg_sum_square / total_sum_square  
    
    return(list(beta    = beta,
                residual          = residual, 
                r_square          = r_square,
                RSS               = residual_sum_squre))
  }
  else { print("ERROR! group not match data!")
  }
}
regfitPermu = function(X, y,group, nperm) {
  #Desp: function to calculate the significance of one fitting
  y               =  as.matrix(y)
  X               =  as.matrix(X)  
  fit_grpreg      =  regfit(X, y, group)
  count_perm      = 0 
  r_square_perm   = seq_len(nperm)
  beta_permu      = as.data.frame(matrix(0, nrow=nperm, ncol = length(fit_grpreg$beta))) 
  while (count_perm < as.numeric(nperm) ) {
    y_permu       =  y + sample(fit_grpreg$residual,size=length(y),replace=TRUE)
    reg_permu     =  regfit(X, y_permu, group)
    r_square_perm[count_perm] = reg_permu$r_square
    beta_permu[count_perm,]   =  reg_permu$beta
    count_perm    =   count_perm  + 1
    if(count_perm > 100 && (count_perm %% 100) == 0) {
      print(paste("Permutation", count_perm))}
  }
  r_square_perm   = r_square_perm[-1]
  pvalue.r2       = max(1/nperm, 
                          length(r_square_perm[which(r_square_perm > fit_grpreg$r_square)]) / nperm)
  bpval = function(x, list){
    xbar = mean(list)
    z <- (xbar-x)/(sd(list)/sqrt(length(list)))
    return(2*pnorm(-abs(z)))
  }
  pvalue          = vapply(seq_along(fit_grpreg$beta), FUN=function(x){bpval(fit_grpreg$beta[x],beta_permu[,x])},0.2)
  names(pvalue)   = names(fit_grpreg$beta)
  
  return(list(RSS      = fit_grpreg$RSS,
              fit_residual = fit_grpreg$residual,
              fit_r2   = c(fit_grpreg$r_square,pvalue = pvalue.r2),
              beta     = cbind(beta = fit_grpreg$beta,pvalue),
              permu_r2 = r_square_perm,
              grpFit   = fit_grpreg ))
  
}
getCandi    = function(fitPermuBeta, pcut){
  x = fitPermuBeta[which(fitPermuBeta[,2] <= pcut),]
  if (dim(x)[1] > 0){
    x = x[x[,1]>0,]
    x = x[order(x[,1],decreasing=T),]
    return(x)
  }else{
    return(c("", ""))
  }
}

#---funce
nperm     = 1000
pvalcut   = 0.01
##---init---

# #----getting command line parameters
args = getArgs()
if(length(args) < 1 || is.null(args)){
  print(paste(usage,example,sep="\n"))
  print(args)
  stop(paste(error,"wrong input parameter!"))
}

setwd(system("pwd",intern=T))
fexp        = args$input
output      = args$output
print(fexp)
print(paste("inputfile",fexp,class(fexp)))
print(paste("outputfile",output))
# output      = jxy(fexp,"_candiReg_", CDT,".txt")


##loading data
rawExp = read.table(fexp,header=T) 
target = as.character(rawExp$gene[1])
regulators = vapply(rawExp$gene[-1],as.character,'a')
samples = colnames(rawExp)[-1]
numReg = length(regulators)
numSmps = length(samples)
regExp = as.data.frame(t(apply(rawExp[,-1],c(1,2),as.numeric)))
rownames(regExp) = samples
colnames(regExp) = c(target, regulators)
regExp = regExp[order(regExp[,1]),]

dataExp = regExp
mycolor = colorRampPalette(c("blue","white","red"))(255)
dataExp = apply(dataExp,2,normalize)
rownames(dataExp) = samples

if (numReg > 2){
  require(adegenet)
  require(amap)
  exp = t(dataExp[,-1])
  expcls = find.clusters(exp,stat='AIC',n.pca=floor(numReg/2),criterion="min",choose.n.clust=FALSE)
  groups = as.character(expcls$grp)
  names(groups) = names(expcls$grp)
}else{
  groups = rep(1, length(regulators))
}


require(grpreg)
print("Doing regression...")
fitpermut = regfitPermu(X=dataExp[,-1],y=dataExp[,1], groups, nperm)

candi = getCandi(fitpermut$beta,pcut=pvalcut)
out = rbind(c("#target", target),
            c("#totalReg",numReg),
            c("#sigReg", nrow(candi)),
            cbind(c("#r2","#r2.pval"),fitpermut$fit_r2),
            c("#regulator\tcoeff","pvalue")
            )

write.table(as.matrix(out),
            output,
            quote=F,col.names=F,sep="\t",row.names = F)
write.table(candi,
            output,
            append=T,
            col.names=F, quote=F,sep="\t")
print(paste("#-Done--",target))

