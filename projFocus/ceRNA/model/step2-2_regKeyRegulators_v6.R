#!/ifs/home/c2b2/ac_lab/jh3283/tools/R/R-3-02/bin/Rscript
#Author: Jing He
#input: < expression matrix (first row as the target, others as regulator) > 
## including fdr p-valu adjustment 

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
##---test
# fexp = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/temp-genelist.failed_Apr11.txt/CDK12_candidateRegs_Apr-11-2014.txt_reg_exp_tumor.temp"
# setwd("/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/temp-genelist.failed_Apr11.txt/")
# output  = "/Volumes/ifs/data/c2b2/ac_lab/jh3283/projFocus/result/03102014/candiReg/temp-genelist.failed_Apr11.txt/CDK12_candidateRegs_Apr-11-2014.local.txt"
##---test

#---funcs
regfit      = function(X, y, group){
  if (setequal(colnames(X), names(group))){
        rss   = NULL
        a     = seq(0.01,1,0.1)
        for (alpha in a){
          grpfit            = grpreg(X,y,group=group,penalty="grLasso",
                                     alpha = alpha,eps=0.05,max.iter=5000)
          grpfit_selected   = select(grpfit,criterion="AIC")
          rss               = c(rss,sum(y - grpfit_selected$beta[1] -
                                          X %*% as.matrix(grpfit_selected$beta[-1])))
        }
        alpha = a[which.min(rss)[1]]
    #fit model
    grpfit              = grpreg(X, y, group=group,penalty="grLasso",
                                  alpha = alpha,eps=0.05,max.iter=5000)
#                                  eps=0.05,max.iter=5000)
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
bpval = function(x, list){
  xbar = mean(list)
  z <- (xbar-x)/(sd(list)/sqrt(length(list)))
  return(2*pnorm(-abs(z)))
}
regfitPermu = function(X, y,group, nperm) {
  #Desp: function to calculate the significance of one fitting
  y               =  as.matrix(y)
  X               =  as.matrix(X)  
  nperm           =  as.numeric(nperm)
  lenY            =   length(y)
  print("Doing first fitting...")
  fit_grpreg      =  regfit(X, y, group)
#   save(fit_grpreg,numReg,target,regulators, file=jxy(output,"permut-temp.rda"))
  count_perm      = 0 
  r_square_perm   = seq_len(nperm)
  beta_permu      = as.data.frame(matrix(0, nrow=nperm, ncol = length(fit_grpreg$beta)))
  print("Start permutation ...")
  pb              = txtProgressBar(min=1, max = nperm, style = 3)
  while (count_perm < nperm ) {
    y_permu                   =  y + sample(fit_grpreg$residual,size=lenY,replace=TRUE)
    reg_permu                 =  regfit(X, y_permu, group)
    r_square_perm[count_perm] = reg_permu$r_square
    beta_permu[count_perm,]   =  reg_permu$beta
    count_perm    =   count_perm  + 1
    setTxtProgressBar(pb,count_perm)
    if((count_perm %% 10) == 0) {
      print(paste("Permutation", count_perm))
    }
  }
  r_square_perm   = r_square_perm[-1]
  pvalue.r2       = max(1/nperm, 
                        length(r_square_perm[which(r_square_perm > fit_grpreg$r_square)]) / nperm)
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
  xname = rownames(fitPermuBeta)
# xname = rownames(fitPermuBeta)[which(fitPermuBeta[,2] <= pcut)]
# x = as.data.frame(fitPermuBeta[which(fitPermuBeta[,2] <= pcut),])
  x = as.data.frame(fitPermuBeta)
  rownames(x) = xname
  if (length(x[,1]) > 1){
    x = x[x[,1]>0,]
    x = x[order(x[,1],decreasing=T),]
    return(x)
  }else if( length(x[,1]) == 1){
    return(x)
  }else {
    return(c("", ""))
  }
}

#---func
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
nperm       = as.integer(args$nperm)
output      = args$output
print(fexp)
print(paste("inputfile",fexp,class(fexp)))
print(paste("outputfile",output))
#---init

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
  groups = as.numeric(as.character(expcls$grp))
  names(groups) = names(expcls$grp)
  if(length(unique(groups)) == 1){
      groups = 1:length(groups)
      names(groups) = regulators
  }
}else if(numReg == 2){
  groups = 1:length(regulators)
  names(groups) = regulators
}else{
  fitlm  = lm(dataExp[,1]~dataExp[,-1])
  sumfit = summary(fitlm)
  save(sumfit,target, regulators, file=jxy(output,".rda"))
  if (sumfit$coefficients[2,1] > 0 ){
    out = rbind(c("#target", target),
                c("#totalReg",1),
                c("#sigReg", 1),
                c("#r2",sumfit$r.squared),
                c("#regulator\tcoeff","pvalue"),
                c(paste(regulators[1], sumfit$coefficients[2,1],sep="\t"),sumfit$coefficients[2,4]) 
    )  
  }else{
    out = rbind(c("#target", target),
                c("#totalReg",1),
                c("#sigReg",0),
                c("#r2",sumfit$r.squared),
                c("#regulator\tcoeff","pvalue")
    )
  }
  write.table(as.matrix(out),
              output,
              quote=F,col.names=F,sep="\t",row.names = F)
  q(save="no")
}

require(grpreg)
print("Doing regression...")
print(table(groups))
fitpermut = regfitPermu(X=dataExp[,-1],y=dataExp[,1], groups, nperm)

save(fitpermut,numReg,target, regulators, file=jxy(output,".rda"))

resBeta = fitpermut$beta
resBeta[,2] = p.adjust(fitpermut$beta[,2],method="fdr")
candi = getCandi(resBeta,pcut=pvalcut)
# candi = resBeta[resBeta[,1]>0,]

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

